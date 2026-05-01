import os
import glob

STRIDE_BIN = os.environ.get("STRIDE_BIN", "stride")
PDB_DIR = os.environ.get("PDB_DIR", "pdbs")
PDBS = [os.path.basename(f).replace(".pdb.gz", "") for f in glob.glob(f"{PDB_DIR}/*.pdb.gz")]

rule all:
    input:
        expand("stride_output/{pdb}.ss.out", pdb=PDBS)

rule unzip_pdb:
    input:
        pdb = f"{PDB_DIR}/{{pdb}}.pdb.gz"
    output:
        unzipped = temp("unzipped_pdbs/{pdb}.pdb")
    shell:
        "zcat {input.pdb} > {output.unzipped}"

rule run_stride:
    input:
        unzipped = "unzipped_pdbs/{pdb}.pdb"
    output:
        ss_out = "stride_output/{pdb}.ss.out"
    shell:
        "{STRIDE_BIN} {input.unzipped}> {output.ss_out}"
