import pandas as pd

import os

# === SETTINGS ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")


# === INPUT / OUTPUT ===
INPUT_FILE = "trf_filtered.tsv"
OUTPUT_FILE = "telomere_candidates.tsv"


# === UTILITY FUNCTIONS ===
def rotations(seq):
    seq = seq.upper()
    return {seq[i:] + seq[:i] for i in range(len(seq))}

def revcomp(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()

# === BUILD ACCEPTED MOTIF SET ===
TEL_MOTIF = "TTAGGG"
TEL_RC = revcomp(TEL_MOTIF)
print(TEL_RC)
all_rotations = rotations(TEL_MOTIF).union(rotations(TEL_RC))
print(f"[INFO] Accepting these telomere motif rotations: {sorted(all_rotations)}")

# === LOAD FILTERED TRF TABLE ===
print(f"[INFO] Loading {INPUT_FILE}...")
df = pd.read_csv(INPUT_FILE, sep="\t")
print(f"[INFO] Loaded {len(df)} entries.")

# === MATCH TEL MOTIFS ===
print("[INFO] Filtering for telomere-like motifs...")
df["motif_upper"] = df["motif"].str.upper()
df_telomere = df[df["motif_upper"].isin(all_rotations)].copy()

print(f"[INFO] Found {len(df_telomere)} telomere-matching entries.")

# === Apply filters ===
min_total_length = 60
min_percent_match = 90

df_filtered = df_telomere[
    (df_telomere["total_length"] >= min_total_length) &
    (df_telomere["percent_matches"] >= min_percent_match)
].copy()

print(f"[INFO] Retained {len(df_filtered)} arrays after filtering.")

print(df_filtered[["chrom", "start", "end", "total_length", "copy_number", "percent_matches", "motif"]].head())

# if you have a .fai file
fai = pd.read_csv("primary_dedup_chr_masked_hp_sealed.fa.fai", sep="\t",usecols=[0,1] ,header=None, names=["chrom","length"])
sizes = dict(zip(fai["chrom"], fai["length"]))


print(sizes)
# define 10kb end zones
df_filtered["near_end"] = df_filtered.apply(
    lambda row: row["start"] <= 10000 or (sizes.get(row["chrom"], 0) - row["end"] <= 10000),
    axis=1
)

print("Arrays near scaffold ends:", df_filtered["near_end"].sum())

# === SAVE ===
df_filtered.drop(columns=["motif_upper"]).to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"[DONE] Saved telomere candidates to {OUTPUT_FILE}")
