from pathlib import Path
import csv
import traceback

import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
import numpy as np

OUT_DIR = Path("/scratch/wvirany/boltz-UQ/output/rnp_full")
GT_DIR  = Path("/scratch/wvirany/ground_truth")
RESULTS = Path("/scratch/wvirany/boltz-UQ/data/results_per_residue.csv")

FIELDNAMES = ["system_id", "chain_id", "res_id", "plddt", "lddt"]


def compute_lddt(pred_cif_path, gt_cif_path, system_id):
    """Compute per-residue lDDT for protein chain(s), aligned to GT receptor.cif.
    
    Handles:
    - Homodimers: Protenix may predict extra chains; we filter pred to only
      the N protein chains specified in the system_id
    - Multi-chain proteins: res_id intersection is done per-chain to avoid
      cross-chain collisions
    """
    protein_chains = [c.split(".")[1] for c in system_id.split("__")[2].split("_")]

    pred_f = pdbx.CIFFile.read(pred_cif_path)
    gt_f = pdbx.CIFFile.read(gt_cif_path)

    pred = pdbx.get_structure(pred_f, model=1, extra_fields=["b_factor"])
    gt = pdbx.get_structure(gt_f, model=1)

    pred_ca = pred[struc.filter_amino_acids(pred) & (pred.atom_name == "CA")]
    gt_ca = gt[struc.filter_amino_acids(gt) & (gt.atom_name == "CA")]

    pred_ca = pred_ca[np.isin(pred_ca.chain_id, protein_chains)]

    gt_chains = list(dict.fromkeys(gt_ca.chain_id))

    all_res_ids, all_plddt, all_lddt, all_chains = [], [], [], []

    for pred_chain, gt_chain in zip(protein_chains, gt_chains):
        p = pred_ca[pred_ca.chain_id == pred_chain]
        g = gt_ca[gt_ca.chain_id == gt_chain]

        common = np.intersect1d(p.res_id, g.res_id)
        if len(common) == 0:
            continue

        p = p[np.isin(p.res_id, common)]
        g = g[np.isin(g.res_id, common)]

        lddt_scores = struc.lddt(g, p, aggregation="residue")

        all_res_ids.append(common)
        all_plddt.append(p.b_factor)
        all_lddt.append(lddt_scores)
        all_chains.append(np.array([pred_chain] * len(common)))

    if not all_res_ids:
        raise ValueError(f"No common residues found between pred and GT for any chain")

    return (
        np.concatenate(all_res_ids),
        np.concatenate(all_plddt),
        np.concatenate(all_lddt),
        np.concatenate(all_chains)
    )


# --- resume logic ---
done = set()
if RESULTS.exists():
    with open(RESULTS) as f:
        reader = csv.DictReader(f)
        for row in reader:
            done.add(row["system_id"])
    print(f"Resuming: {len(done)} systems already done")

write_header = not RESULTS.exists()

system_dirs = sorted(OUT_DIR.iterdir())

with open(RESULTS, "a", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    for i, system_dir in enumerate(system_dirs):
        system_id = system_dir.name
        if system_id in done:
            continue

        pred_cif = system_dir / "seed_101" / "predictions" / f"{system_id}_sample_0.cif"
        gt_cif   = GT_DIR / system_id / "receptor.cif"

        if not pred_cif.exists() or not gt_cif.exists():
            print(f"SKIP {system_id}: missing file")
            continue

        try:
            common_res, plddt, lddt, chains = compute_lddt(pred_cif, gt_cif, system_id)
            for j in range(len(common_res)):
                writer.writerow({
                    "system_id":  system_id,
                    "chain_id": chains[j],
                    "res_id":     int(common_res[j]),
                    "plddt":      float(plddt[j]),
                    "lddt":       float(lddt[j]),
                })
            f.flush()
            if i % 100 == 0:
                print(f"[{i}/{len(system_dirs)}] OK {system_id}: {len(common_res)} residues")
        except Exception as e:
            print(f"ERROR {system_id}: {e}")
            traceback.print_exc()