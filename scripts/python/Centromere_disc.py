import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import subprocess
import networkx as nx

# === SETTINGS ===
os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === INPUT / OUTPUT FILES ===
INPUT_FILE = "trf_filtered.tsv"
FILTERED_TSV = "centromere_candidates.tsv"
FASTA_FILE = "centromere_candidates.fasta"

# === 1. Load TRF-filtered table ===
df = pd.read_csv(INPUT_FILE, sep="\t")

# === 2. Apply centromere candidate filter ===
df_centromere = df[
    (df["total_length"] >= 10_000) &
    (df["period"] >= 150) &
    (df["copy_number"] >= 10) &
    (df["percent_matches"] >= 75)
].copy()

print(f"[INFO] Found {len(df_centromere)} candidate centromeric arrays.")

# === 3. Save filtered table ===
df_centromere.to_csv(FILTERED_TSV, sep="\t", index=False)

# === 4. Write motif sequences to FASTA ===
seq_records = [
    SeqRecord(Seq(row["motif"].strip()),
              id=f"{row['chrom']}:{row['start']}-{row['end']}",
              description="")
    for _, row in df_centromere.iterrows()
]

SeqIO.write(seq_records, FASTA_FILE, "fasta")

print(f"[INFO] Wrote FASTA file with {len(seq_records)} consensus repeat units → {FASTA_FILE}")
print(f"[INFO] Wrote filtered table → {FILTERED_TSV}")


BLAST_DB_PREFIX = "centromere_db"
BLAST_RAW_OUTPUT = "centromere_blast.tsv"
FILTERED_OUTPUT = "centromere_filtered_similarity.tsv"

MAKEBLASTDB_CMD = "/Users/lizicheng/opt/anaconda3/envs/biogpt_dl/bin/makeblastdb"
BLASTN_CMD = "/Users/lizicheng/opt/anaconda3/envs/biogpt_dl/bin/blastn"

# === 1. Make BLAST database ===
print("[INFO] Creating BLAST database...")
subprocess.run([
    MAKEBLASTDB_CMD,
    "-in", FASTA_FILE,
    "-dbtype", "nucl",
    "-out", BLAST_DB_PREFIX
], check=True)

# === 2. Run all-vs-all BLASTN ===
print("[INFO] Running all-vs-all BLAST...")
subprocess.run([
    BLASTN_CMD,
    "-query", FASTA_FILE,
    "-db", BLAST_DB_PREFIX,
    "-outfmt", "6 qseqid sseqid pident length qlen slen",
    "-task", "blastn",
    "-evalue", "1e-5",
    "-out", BLAST_RAW_OUTPUT
], check=True)

# === 3. Filter BLAST results ===
print("[INFO] Filtering BLAST hits...")
cols = ["query", "subject", "pident", "length", "qlen", "slen"]
blast_df = pd.read_csv(BLAST_RAW_OUTPUT, sep="\t", names=cols)

# Remove self-matches
blast_df = blast_df[blast_df["query"] != blast_df["subject"]].copy()

# Compute coverage over the shorter sequence
blast_df["min_len"] = blast_df[["qlen", "slen"]].min(axis=1)
blast_df["coverage"] = blast_df["length"] / blast_df["min_len"]

# Apply thresholds
filtered_df = blast_df[
    (blast_df["pident"] >= 80) &
    (blast_df["coverage"] >= 0.7)
].copy()

# Keep only top-scoring hit per query–subject pair
filtered_df = (
    filtered_df
    .sort_values(by=["query", "subject", "pident", "length"], ascending=[True, True, False, False])
    .drop_duplicates(subset=["query", "subject"])
)

# === 4. Save filtered table ===
filtered_df.to_csv(FILTERED_OUTPUT, sep="\t", index=False)
print(f"[INFO] Filtered similarity table written to: {FILTERED_OUTPUT}")
print(f"[INFO] Retained {len(filtered_df)} high-confidence hits.")



FILTERED_BLAST_FILE = f"centromere_filtered_similarity.tsv"
CLUSTERS_OUTPUT_FILE = f"centromere_clusters.tsv"

# === 1. Load filtered BLAST similarity table ===
print("[INFO] Loading filtered similarity file...")
df = pd.read_csv(FILTERED_BLAST_FILE, sep="\t")
print(f"[INFO] Loaded {len(df)} pairwise hits.")

# === 2. Build undirected graph from query-subject pairs ===
print("[INFO] Building similarity graph...")
G = nx.Graph()
G.add_edges_from(zip(df["query"], df["subject"]))

# === 3. Find connected components (clusters) ===
print("[INFO] Finding connected components...")
clusters = list(nx.connected_components(G))
print(f"[INFO] Found {len(clusters)} clusters.")

# === 4. Assign cluster IDs ===
cluster_table = []
for cluster_id, members in enumerate(clusters, 1):
    for seq_id in members:
        cluster_table.append((seq_id, cluster_id))

# === 5. Save cluster assignments ===
cluster_df = pd.DataFrame(cluster_table, columns=["seq_id", "cluster_id"])
cluster_df.to_csv(CLUSTERS_OUTPUT_FILE, sep="\t", index=False)
print(f"[INFO] Cluster table written to {CLUSTERS_OUTPUT_FILE}")