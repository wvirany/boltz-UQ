import json
import os
import string

inputs_json = "/scratch/wvirany/inputs.json"
msa_dir = "/scratch/wvirany/msa_files"
output_dir = "/scratch/wvirany/boltz-UQ/inputs/rnp_full"

os.makedirs(output_dir, exist_ok=True)

data = json.load(open(inputs_json))
chain_ids = list(string.ascii_uppercase)

skipped = []
for system_id, system_data in data.items():
    msa_system_dir = os.path.join(msa_dir, system_id.lower())
    if not os.path.exists(msa_system_dir):
        skipped.append(system_id)
        continue

    sequences_list = []
    for idx, (chain, seq) in enumerate(system_data["sequences"].items()):
        chain_id = chain_ids[idx]
        sequences_list.append({
            "proteinChain": {
                "sequence": seq,
                "count": 1,
                "msa": {
                    "precomputed_msa_dir": os.path.join(msa_system_dir, chain_id),
                    "pairing_db": "uniprot"
                }
            }
        })

    for smiles in system_data["smiles"]:
        sequences_list.append({
            "ligand": {
                "ligand": smiles,
                "count": 1
            }
        })

    output = [{"name": system_id, "sequences": sequences_list}]
    with open(os.path.join(output_dir, f"{system_id}.json"), "w") as f:
        json.dump(output, f, indent=2)

print(f"Written: {len(data) - len(skipped)}")
print(f"Skipped: {len(skipped)}")
if skipped:
    print("Skipped IDs:", skipped[:5], "...")