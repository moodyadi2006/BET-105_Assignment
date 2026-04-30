import sys
import pandas as pd
import re
import os

stride_output_file = sys.argv[1]
output_file = sys.argv[2]
target_aa = sys.argv[3]

pdb_id = os.path.basename(stride_output_file).split(".")[0]

three_to_one_letter = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

residue_records = []

with open(stride_output_file) as f:
    for line in f:
        if line.startswith("ASG"):
            parts = line.split()

            if len(parts) < 10:
                continue

            resnum_match = re.match(r"\d+", parts[3])
            if not resnum_match:
                continue

            try:
                residue_records.append({
                    "res3": parts[1],
                    "chain": parts[2],
                    "resnum": int(resnum_match.group()),
                    "idx": int(parts[4]),
                    "ss_code": parts[5],
                    "ss_name": parts[6],
                    "phi": float(parts[7]),
                    "psi": float(parts[8]),
                    "area": float(parts[9]),
                })
            except:
                continue

residue_df = pd.DataFrame(residue_records)

tripeptide_rows = []

for i in range(1, len(residue_df) - 1):
    center_res = residue_df.iloc[i]

    if center_res["res3"] != target_aa:
        continue

    prev_res = residue_df.iloc[i - 1]
    next_res = residue_df.iloc[i + 1]

    tripeptide_seq = (
        three_to_one_letter.get(prev_res["res3"], "X") +
        three_to_one_letter.get(center_res["res3"], "X") +
        three_to_one_letter.get(next_res["res3"], "X")
    )
    ss_triplet = prev_res["ss_code"] + center_res["ss_code"] + next_res["ss_code"]

    position_string = (
        f"{center_res['chain']}:{prev_res['resnum']},"
        f"{center_res['chain']}:{center_res['resnum']},"
        f"{center_res['chain']}:{next_res['resnum']}"
    )

    for res in [prev_res, center_res, next_res]:
        tripeptide_rows.append([
            res["res3"],
            res["chain"],
            res["resnum"],
            res["idx"],
            res["ss_code"],
            res["ss_name"],
            res["phi"],
            res["psi"],
            res["area"],
            three_to_one_letter.get(res["res3"], "X"),
            tripeptide_seq,
            ss_triplet,
            pdb_id,
            position_string
        ])

output_df = pd.DataFrame(tripeptide_rows)
output_df.to_csv(output_file, sep="\t", index=False, header=False)
