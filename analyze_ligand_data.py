from pathlib import Path
import csv
import json
import traceback

import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

OUT_DIR    = Path("/scratch/wvirany/boltz-UQ/output/rnp_full")
GT_DIR     = Path("/scratch/wvirany/ground_truth")
INPUT_DIR  = Path("/scratch/wvirany/boltz-UQ/inputs/rnp_full")
RESULTS    = Path("/scratch/wvirany/boltz-UQ/results_ligand_rmsd.csv")

FIELDNAMES = ["system_id", "ligand_chain_id", "n_atoms", "plddt_mean", "rmsd"]


def compute_ligand_rmsd(pred_cif_path, gt_dir, input_json_path, system_id):
    ligand_chains = system_id.split("__")[3].split("_")

    f = pdbx.CIFFile.read(pred_cif_path)
    s = pdbx.get_structure(f, model=1, extra_fields=["b_factor"])

    all_chains = list(dict.fromkeys(s.chain_id))
    protein_chains = [c.split(".")[1] for c in system_id.split("__")[2].split("_")]
    pred_ligand_chains = [c for c in all_chains if c not in protein_chains]

    with open(input_json_path) as f:
        inp = json.load(f)
    ligand_entries = [e["ligand"]["ligand"] for e in inp[0]["sequences"] if "ligand" in e]

    rows = []
    for pred_chain, gt_chain_id, smiles in zip(pred_ligand_chains, ligand_chains, ligand_entries):
        atoms = s[s.chain_id == pred_chain]

        if len(atoms) <= 1:
            continue

        pred_coords   = atoms.coord
        pred_elements = [e.capitalize() for e in atoms.element]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not parse SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mol = Chem.RemoveHs(mol)

        mol_elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
        if pred_elements != mol_elements:
            raise ValueError(f"Atom order mismatch for {gt_chain_id}: "
                             f"{pred_elements} vs {mol_elements}")

        conf = mol.GetConformer()
        for j in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(j, pred_coords[j].tolist())

        gt_sdf = gt_dir / system_id / "ligand_files" / f"{gt_chain_id}.sdf"
        gt_mol = Chem.MolFromMolFile(str(gt_sdf), removeHs=True)
        if gt_mol is None:
            raise ValueError(f"Could not load GT SDF: {gt_sdf}")

        rmsd = rdMolAlign.GetBestRMS(gt_mol, mol)
        rows.append({
            "system_id":       system_id,
            "ligand_chain_id": gt_chain_id,
            "n_atoms":         len(atoms),
            "plddt_mean":      float(atoms.b_factor.mean()),
            "rmsd":            float(rmsd),
        })

    return rows


# --- resume logic ---
done = set()
if RESULTS.exists():
    with open(RESULTS) as f:
        reader = csv.DictReader(f)
        for row in reader:
            done.add(row["system_id"])
    print(f"Resuming: {len(done)} systems already done")

write_header = not RESULTS.exists()
system_dirs  = sorted(OUT_DIR.iterdir())

with open(RESULTS, "a", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    for i, system_dir in enumerate(system_dirs):
        system_id = system_dir.name
        if system_id in done:
            continue

        pred_cif   = system_dir / "seed_101" / "predictions" / f"{system_id}_sample_0.cif"
        input_json = INPUT_DIR / f"{system_id}-update-msa.json"

        if not pred_cif.exists():
            print(f"SKIP {system_id}: no predicted CIF")
            continue
        if not input_json.exists():
            print(f"SKIP {system_id}: no input JSON")
            continue

        try:
            rows = compute_ligand_rmsd(pred_cif, GT_DIR, input_json, system_id)
            for row in rows:
                writer.writerow(row)
            f.flush()
            if i % 100 == 0:
                print(f"[{i}/{len(system_dirs)}] OK {system_id}: {len(rows)} ligands")
        except Exception as e:
            print(f"ERROR {system_id}: {e}")
            traceback.print_exc()