import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import os
import networkx as nx

# === Paths ===
dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/"
input_tsv = os.path.join(dir, "centromere_candidates.tsv")
output_dir = os.path.join(dir, "blast")
os.makedirs(output_dir, exist_ok=True)

fasta_file = os.path.join(output_dir, "centromeres.fasta")
blast_db_prefix = os.path.join(output_dir, "centromere_db")
blast_output_file = os.path.join(output_dir, "centromere_blast.tsv")
filtered_output = os.path.join(output_dir, "filtered_similarity.tsv")

# === 1. Load candidate centromere sequences ===
df = pd.read_csv(input_tsv, sep="\t")
seq_records = [
    SeqRecord(Seq(row["repeat_seq"]), id=f"{row['chrom']}:{row['start']}-{row['end']}", description="")
    for _, row in df.iterrows()
]

# === 2. Write sequences to FASTA ===
SeqIO.write(seq_records, fasta_file, "fasta")

# === 3. Make BLAST database ===
subprocess.run([
    "/Users/lizicheng/opt/anaconda3/envs/biogpt_dl/bin/makeblastdb",
    "-in", fasta_file, "-dbtype", "nucl", "-out", blast_db_prefix
], check=True)

# === 4. All-vs-all BLASTN ===
subprocess.run([
    "/Users/lizicheng/opt/anaconda3/envs/biogpt_dl/bin/blastn",
    "-query", fasta_file, "-db", blast_db_prefix,
    "-outfmt", "6 qseqid sseqid pident length qlen slen",
    "-task", "blastn", "-evalue", "1e-5", "-out", blast_output_file
], check=True)

# === 5. Filter BLAST results based on similarity and coverage ===
cols = ["query", "subject", "pident", "length", "qlen", "slen"]
blast_df = pd.read_csv(blast_output_file, sep="\t", names=cols)

# Remove self matches
blast_df = blast_df[blast_df["query"] != blast_df["subject"]].copy()

# Compute coverage of the alignment over the shorter sequence
blast_df["min_len"] = blast_df[["qlen", "slen"]].min(axis=1)
blast_df["coverage"] = blast_df["length"] / blast_df["min_len"]

# Apply thresholds
filtered = blast_df[
    (blast_df["pident"] >= 75) &
    (blast_df["coverage"] >= 0.50)
].copy()

# Keep only top-scoring hit per query–subject pair
filtered = filtered.sort_values(
    by=["query", "subject", "pident", "length"],
    ascending=[True, True, False, False]
).drop_duplicates(subset=["query", "subject"])

# === 6. Save and report ===
filtered.to_csv(filtered_output, sep="\t", index=False)
print(f"Done. Retained {len(filtered)} query–subject pairs with ≥75% identity and ≥90% coverage.")


# === Load filtered BLAST similarity table ===
filtered_file = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/blast/filtered_similarity.tsv"
df = pd.read_csv(filtered_file, sep="\t")

# === Build undirected graph from query-subject pairs ===
G = nx.Graph()
G.add_edges_from(zip(df["query"], df["subject"]))

# === Find connected components (clusters) ===
clusters = list(nx.connected_components(G))

# === Output clusters ===
cluster_table = []
for cluster_id, members in enumerate(clusters, 1):
    for seq_id in members:
        cluster_table.append((seq_id, cluster_id))

# === Save cluster membership table ===
cluster_df = pd.DataFrame(cluster_table, columns=["seq_id", "cluster_id"])
out_path = filtered_file.replace("filtered_similarity.tsv", "centromere_clusters.tsv")
cluster_df.to_csv(out_path, sep="\t", index=False)

print(f"Identified {len(clusters)} clusters.")
print(cluster_df)

# # === Define Cluster 2 sequences (replace with your actual set) ===
# cluster2 = {
# "HiC_scaffold_1:35332756-35343344",
# "HiC_scaffold_1:35253705-35280255",
# "HiC_scaffold_8:29605285-29639983",
# "HiC_scaffold_9:20247310-20270222",
# "HiC_scaffold_1:35162571-35213205",
# "HiC_scaffold_9:20171149-20199590",
# "HiC_scaffold_6:22463395-22501347",
# "HiC_scaffold_10:16428140-16477937",
# "HiC_scaffold_8:29680790-29693725",
# "HiC_scaffold_1:35284247-35295986",
# "HiC_scaffold_3:31012437-31062449",
# "HiC_scaffold_6:22422658-22441715",
# "HiC_scaffold_10:16309140-16389363"
# }

# # === Filter all query-subject pairs where BOTH are in Cluster 2 ===
# cluster2_alignments = filtered[
#     (blast_df["query"].isin(cluster2)) & 
#     (blast_df["subject"].isin(cluster2))
# ]

# # === Print or save the result ===
# print(f"Cluster 2 has {len(cluster2_alignments)} alignment pairs:\n")
# print(cluster2_alignments)

# # Optional: Save to file
# cluster2_alignments.to_csv(
#     "/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/stlt_analysis/blast/cluster2_alignments.tsv",
#     sep="\t", index=False
# )
