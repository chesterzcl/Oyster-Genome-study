import pandas as pd
import os

# === SETTINGS ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")


# === PARAMETERS ===
INPUT_FILE = "primary_dedup_chr_masked_hp_sealed.fa.2.7.7.80.10.50.2000.dat"
OUTPUT_FILE = "trf_filtered.tsv"

# === GENERIC QUALITY FILTER THRESHOLDS ===
MIN_PERIOD = 2
MIN_COPY_NUMBER = 2
MIN_CONSENSUS_SIZE = 3
MIN_PERCENT_MATCHES = 70

# === PARSE TRF FILE ===
print("[INFO] Parsing TRF .dat file...")
records = []
current_seq = None

with open(INPUT_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith("Sequence:"):
            current_seq = line.split(":", 1)[1].strip()
            continue
        if line.startswith("Parameters:"):
            continue

        parts = line.split()
        if len(parts) < 12:
            continue

        try:
            start = int(parts[0])
            end = int(parts[1])
            period = int(parts[2])
            copy_number = float(parts[3])
            consensus_size = int(parts[4])
            percent_matches = int(parts[5])
            total_length = end - start + 1
            motif = parts[13] if len(parts) > 13 else "."
        except Exception as e:
            print(f"[WARN] Skipping line: {line}\nReason: {e}")
            continue

        records.append({
            "chrom": current_seq,
            "start": start,
            "end": end,
            "total_length": total_length,
            "period": period,
            "copy_number": copy_number,
            "consensus_size": consensus_size,
            "percent_matches": percent_matches,
            "motif": motif
        })

df = pd.DataFrame(records)
print(f"[INFO] Parsed {len(df)} entries.")

# === APPLY QUALITY FILTER ===
print("[INFO] Applying quality filters...")
df_filtered = df[
    (df["period"] >= MIN_PERIOD) &
    (df["copy_number"] >= MIN_COPY_NUMBER) &
    (df["consensus_size"] >= MIN_CONSENSUS_SIZE) &
    (df["percent_matches"] >= MIN_PERCENT_MATCHES)
].copy()

print(f"[INFO] Remaining after filters: {len(df_filtered)} entries.")

# === REMOVE EXACT DUPLICATES ===
print("[INFO] Removing exact duplicates...")
df_filtered = (
    df_filtered
    .sort_values(["percent_matches", "copy_number"], ascending=[False, False])
    .drop_duplicates(subset=["chrom", "start", "end"])
    .reset_index(drop=True)
)

print(f"[INFO] Remaining after deduplication: {len(df_filtered)} entries.")

# === SAVE RESULT ===
df_filtered.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"[DONE] Filtered TRF table saved to {OUTPUT_FILE}")


