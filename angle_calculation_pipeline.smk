configfile: "config.yaml"

import glob
import os

TARGET_AA = config["target_aa"]

PDB_IDS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("stride_output/*.ss.out")]

rule all:
    input:
        "results/angles.tsv"


rule calculate_angles:
    input:
        expand("tripeptide_contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDB_IDS, aa=[TARGET_AA])
    output:
        "results/angles.tsv"
    params:
        aa=TARGET_AA,
        pdb_dir=os.environ.get("PDB_DIR", "pdbs")
    shell:
        """
        python scripts/calculate_sidechain_angles.py tripeptide_contexts {output} {params.aa} {params.pdb_dir}
        """
