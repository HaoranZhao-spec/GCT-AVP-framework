import pandas as pd
import numpy as np
import re
import argparse

def calc_biophysical_features(seq):
    seq = seq.upper().strip()
    motif_score = 0.0
    if seq.startswith("FF") or seq.startswith("FL"):
        motif_score += 1.5
    if "KLLL" in seq or "AKIL" in seq or "AKVL" in seq:
        motif_score += 1.2
    if seq.endswith("FFS") or seq.endswith("LFS") or seq.endswith("LSS"):
        motif_score += 1.0
    hydro_ratio = sum(seq.count(aa) for aa in "FLIVMA") / len(seq)
    pos_ratio = sum(seq.count(aa) for aa in "KR") / len(seq)
    return motif_score + hydro_ratio * 2.0 + pos_ratio * 1.5

def main(input_file, output_file):
    df = pd.read_csv(input_file)
    seq_col = [c for c in df.columns if 'seq' in c.lower() or 'sequence' in c.lower()][0]
    sequences = df[seq_col].astype(str).tolist()
    bio_scores = np.array([calc_biophysical_features(s) for s in sequences])
    raw_pEC50 = df["Predicted_pEC50"].values
    norm_pEC50 = (raw_pEC50 - raw_pEC50.min()) / (raw_pEC50.max() - raw_pEC50.min() + 1e-8)
    norm_bio = (bio_scores - bio_scores.min()) / (bio_scores.max() - bio_scores.min() + 1e-8)
    df["Composite_Score"] = 0.35 * norm_pEC50 + 0.65 * norm_bio
    df["Predicted_EC50_uM"] = 10 ** (6.0 - df["Predicted_pEC50"])
    df_sorted = df.sort_values(by="Composite_Score", ascending=False).reset_index(drop=True)
    df_sorted.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input csv with Predicted_pEC50 and sequence")
    parser.add_argument("--output", required=True, help="Output composite score csv")
    args = parser.parse_args()
    main(args.input, args.output)