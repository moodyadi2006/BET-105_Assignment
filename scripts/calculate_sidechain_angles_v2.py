import os
import sys
import gzip
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm

# -----------------------------
# INPUTS
# -----------------------------
tripeptide_contexts_dir = sys.argv[1]
output_file = sys.argv[2]
target_aa = sys.argv[3]

pdb_structures_dir = "pdb_structures"


residue_size_map = {
    "K": "Large",  "R": "Bulky",        "D": "Intermediate", "E": "Large",
    "G": "Tiny",   "A": "Tiny",         "V": "Small",        "L": "Intermediate",
    "I": "Intermediate", "M": "Large",  "P": "Small",        "S": "Small",
    "T": "Small",  "N": "Intermediate", "Q": "Large",        "C": "Small",
    "F": "Bulky",  "Y": "Bulky",        "W": "Bulky",        "H": "Large"
}


def signed_angle_3d(v1, v2, axis):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    axis = axis / np.linalg.norm(axis)

    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)

    angle = np.arctan2(np.dot(cross, axis), dot)
    return np.degrees(angle)


def load_structure(pdb_id):
    file_path = os.path.join(pdb_structures_dir, f"{pdb_id}.pdb.gz")
    if not os.path.exists(file_path):
        return None
    parser = PDBParser(QUIET=True)
    with gzip.open(file_path, "rt") as f:
        structure = parser.get_structure(pdb_id, f)
    return structure


def get_residue(structure, chain_id, resnum):
    try:
        for model in structure:
            chain = model[chain_id]
            for res in chain:
                if res.id[1] == resnum:
                    return res
    except:
        return None
    return None


def get_alpha_carbon(res):
    if res and "CA" in res:
        return res["CA"].get_coord()
    return None


def get_sidechain_centroid(res):
    coords = []
    for atom in res:
        if atom.get_name() not in ["N", "CA", "C", "O"]:
            coords.append(atom.get_coord())
    if not coords:
        return None
    return np.mean(coords, axis=0)


angle_results = []

context_files = [f for f in os.listdir(tripeptide_contexts_dir) if f.endswith(".tsv")]

for context_file in tqdm(context_files, desc="Processing contexts"):

    context_path = os.path.join(tripeptide_contexts_dir, context_file)

    if os.path.getsize(context_path) == 0:
        continue

    try:
        context_df = pd.read_csv(context_path, sep="\t", header=None)
    except:
        continue

    if context_df.empty:
        continue

    pdb_id = context_df.iloc[0, 12]
    pdb_structure = load_structure(pdb_id)

    if pdb_structure is None:
        continue

    # process tripeptides (3 rows per context)
    for i in range(0, len(context_df), 3):

        try:
            prev_row = context_df.iloc[i]
            center_row = context_df.iloc[i + 1]
            next_row = context_df.iloc[i + 2]
        except:
            continue

        if center_row[0] != target_aa:
            continue

        if center_row[11] != "HHH":
            continue

        chain_id = center_row[1]

        prev_resnum = int(prev_row[2])
        center_resnum = int(center_row[2])
        next_resnum = int(next_row[2])

        res_prev = get_residue(pdb_structure, chain_id, prev_resnum)
        res_cent = get_residue(pdb_structure, chain_id, center_resnum)
        res_next = get_residue(pdb_structure, chain_id, next_resnum)

        if not res_prev or not res_cent or not res_next:
            continue

        ca_prev = get_alpha_carbon(res_prev)
        ca_cent = get_alpha_carbon(res_cent)
        ca_next = get_alpha_carbon(res_next)

        centroid_prev = get_sidechain_centroid(res_prev)
        centroid_cent = get_sidechain_centroid(res_cent)

        if any(x is None for x in [ca_prev, ca_cent, ca_next, centroid_prev, centroid_cent]):
            continue

        sidechain_vec_prev = centroid_prev - ca_prev
        sidechain_vec_cent = centroid_cent - ca_cent

        # alternative axis: direct vector between consecutive alpha carbons
        helix_axis = ca_cent - ca_prev

        if np.linalg.norm(helix_axis) == 0:
            continue

        try:
            sidechain_angle = signed_angle_3d(sidechain_vec_prev, sidechain_vec_cent, helix_axis)
        except:
            continue

        prev_aa_one_letter = prev_row[9]
        size_class = residue_size_map.get(prev_aa_one_letter, "Unknown")

        if size_class == "Unknown":
            continue

        angle_results.append([
            pdb_id,
            prev_aa_one_letter,
            size_class,
            sidechain_angle
        ])

results_df = pd.DataFrame(angle_results, columns=[
    "pdb", "left_aa", "size_class", "angle"
])

os.makedirs(os.path.dirname(output_file), exist_ok=True)
results_df.to_csv(output_file, sep="\t", index=False)

print("Saved:", output_file)
print("Total angles:", len(results_df))
