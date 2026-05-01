#!/bin/bash

CORES=4


export PDB_DIR="${1:-pdbs}"
export STRIDE_BIN="${2:-/usr/local/bin/stride}"


echo " PDB to STRIDE"
snakemake -s secondary_structure_pipeline.smk --cores $CORES --keep-going

echo "STRIDE to Contexts"
snakemake -s context_extraction_pipeline.smk --cores $CORES --keep-going

echo "Contexts to Angle list"
snakemake -s angle_calculation_pipeline.smk --cores $CORES --keep-going

echo "Plot"
Rscript scripts/plot_angle_distribution.R

echo "Done. Output: results/angle_plot.png"
