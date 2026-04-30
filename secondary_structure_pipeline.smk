rule unzip_pdb:
	input:
		pdb = "pdb_structures/{pdb}.pdb.gz"
	output:
		unzipped = temp("unzipped_structures/{pdb}.pdb")
	shell:
		" zcat {input.pdb} > {output.unzipped} "
rule run_stride:
	input:
		unzipped = "unzipped_structures/{pdb}.pdb"
	output:
		ss_out = "stride_output/{pdb}.ss.out"
	params:
	shell:
		"/mnt/c/Users/ADITYA/stride/stride {input.unzipped} > {output.ss_out}"
