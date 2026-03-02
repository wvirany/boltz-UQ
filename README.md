a small project for studying uncertainty quantification in protein structure and binding affinity prediction models

Addtional useful documentation:
- [Main Protenix documentation](https://github.com/bytedance/Protenix/tree/main/docs)
- [DRAC documentation on JupyterLab](https://docs.alliancecan.ca/wiki/JupyterLab#The_JupyterLab_interface)

# Environment setup

The following includes instructions for setting up a Python environment to run Protenix on a DRAC node (tested on Nibi).

**Load modules:**
```bash
module load python/3.11 rdkit/2024.09.6 cuda/12.6
```

**Virtual environment:** Setup (or load) virtual environment in root directory
```bash
python -m venv venv    # If not already created
source venv/bin/activate
pip install --no-index --upgrade pip
```

**Install packages:**
```bash
pip install --no-index torch numpy matplotlib biotite biopython triton optree
pip install --no-index cuequivariance-ops-torch-cu12
pip install fair-esm click ml_collections

pip install --no-deps protenix
```

# Inference

```bash
protenix predict --input {input.json} --out_dir ./output
```

**Args:**
- `--input`: Path to the input JSON file specifying the structure(s) to predict
- `--out_dir`: Output directory for predictions
- `--model_name`: Model variant to use. Some likely options:
	- `protenix_tiny_default_v0.5.0`
	- `protenix_mini_default_v0.5.0`
	- `protenix_base_default_v1.0.0`
- `--use_msa`: Whether to use Multiple Sequence Alignments. Requires external databases if enabled. Set to `false` for quick runs
- `--use_default_params`: When `true`, automatically sets `--cycle` and `--step` based on the model. Set to `false` to override manually
- `--step`: Number of diffusion steps. Default for tiny is 1; more steps results in higher quality predictions but slower
- `--cycle`: Number of recycle iterations through the model. More cycles results in higher quality predictions but slower
- `--sample`: Number of structure samples to generate per input
- `--trimul_kernel` / `--triatt_kernel`: Kernel implementations for triangular multiplicative and attention modules. Use `torch` for compatibility (e.g. on MIG slices), `cuequivariance` for speed on full H100 nodes. See `kernels.md` in main docs for more info.

## Inference example

An example of predicting the structure of [1PIN](https://www.rcsb.org/structure/1PIN):

1. **Download the FASTA file containing the sequence:**
```fasta
>1PIN_1|Chain A|PEPTIDYL-PROLYL CIS-TRANS ISOMERASE|Homo sapiens (9606)
MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQSRRPSSWRQEKITRTKEEALELINGYIQKIKSGEEDFESLASQFSDCSSAKARGDLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHIILRTE
```

2. **Specify the input JSON:**
```
[{
    "name": "test_tiny",
    "sequences": [
        {
            "proteinChain": {
                "sequence": "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQSRRPSSWRQEKITRTKEEALELINGYIQKIKSGEEDFESLASQFSDCSSAKARGDLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHIILRTE",
                "count": 1
            }
        }
    ]
}]
```

This is organized as a list of jobs: `[{job1}, {job2}, ...]`, where each job is a dictionary specifying the job name, sequences of each structure (multiple structures for complexes), and potentially other metadata.

3. **Make a prediction:**
```bash
protenix pred --input 1PIN.json --out_dir ./output --model_name protenix_tiny_default_v0.5.0 --use_msa false --use_default_params false --step 1 --cycle 1 --sample 1 --trimul_kernel torch --triatt_kernel torch
```

Note I used very conservative parameters here for a quick run.

4. **Analyze results:**
```
├── 1PIN.fasta
├── 1PIN.json
└── output
    └── test_tiny
        └── seed_101
            └── predictions
                ├── test_tiny_sample_0.cif
                └── test_tiny_summary_confidence_sample_0.json
```

The confidence prediction JSON file looks like
```
{
    "plddt": 37.274845123291016,
    "gpde": 4.546759605407715,
    "ptm": 0.2619346082210541,
    "iptm": 0.0,
    ...
}
```

# Data processing

## Protein structure validation set

The validation set consists of 100 single-chain protein structures deposited in the PDB after Protenix's training cutoff (2021-09-30), ensuring an unseen validation set.


Candidate structures were retrieved programmatically via the RCSB PDB search API with the following filters:
- Deposition date after `2021-09-30`
- Experimental method: X-ray diffraction only
- Resolution ≤ 2.5 Å
- Single polymer chain instance (strict monomer)
- Polymer type: protein only (no DNA/RNA)
- No small molecule ligands (nonpolymer entity count = 0)
- Sequence length: 50–500 residues

 This yielded 1311 candidate structures. To ensure sequence diversity, candidates were mapped to RCSB's pre-computed 30% sequence identity clusters. One representative was selected at random per cluster, yielding 563 length-filtered representatives, from which 100 were randomly sampled.
 
 **Files:**
 - `data/protein_validation_set.fasta`: sequences for all 100 proteins