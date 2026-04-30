import glob

configfile: "config.yaml"

TARGET_AA = config["target_aa"]


PDB_IDS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("stride_output/*.ss.out")]


rule all:
	input:
		expand("tripeptide_contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDB_IDS, aa=[TARGET_AA])

rule extract_context:
	input:
		ss = "stride_output/{pdb}.ss.out"
	output:
		tsv = "tripeptide_contexts/context_for_{aa}_in_{pdb}.tsv"
	params:
		aa = "{aa}"
	shell:
		"""
		python scripts/extract_tripeptide_context.py {input.ss} {output.tsv} {params.aa}
		"""
