from Bio import pairwise2
from Bio.Seq import Seq
import numpy as np
import pandas as pd

# === 1. Your sequence ===
seq = ("CATTTCCAGAAACAGATCATATTGGCTCAGTGGATAGCGCACTGGACTTGAGATCCAAAGGACGCGAGTTCAAGTCTTACAGCAGCAGATTTTTTTTTCGGAGGGGGGTCTCGTTTTTTTAACATTTATGTATTTAAAAATGTTTGTTCATTTATTTGTACTTTTTTCTTTACACTTGAATTTTTTTTCTATTTTGATGGAAAGCAAATTGTAATTCCAAAATGATGATTAAGACGAAATGCCTCAAACAGGCTAGGATTAAACAAAACGTTTAAAAACGCTTCAAGGGGGAACTTCAGGGGGAGATAACTCACTGCGATATTTGACCAGAATACTAGAGAAAGTAACACTTACATCGTATTGGCCGAGTGAATAACACACTGGATTTGAAAACCAATGGTCGCAATTTCGATTCTTACAGCAGAAGATATTTTTTCTTTGGGCGGGGGTATCGTTTTTTGTCATCATGTATGTGTTTAAAAATGTTTGTTCATTTCTTATTCCTTTTTTCTTAACGTTTGAAAAATTCTTTCCAACTTTGACGGAAAAGGAATTATAATGCCCCAAAAAATTCTTTGAGACGAAATGCCTATAACAAGATTCTAATTTAAAAAAAAGATTAGCCAAATTCAATTTCTCTTTTTGAAAAAGAATGCGTCGGATCAGTTATTTCTCTACACTCACACCAATATAGCAATTTTTTCGAACTGTATATATATGCAGGTGTTTAACATTTATCATTTTCTTGAATTTGAAAACTTGCATATTGTGTATTAATTAACATACGTAATCAAATCGTTTCAGGGGGAGACAACTTCAGGGGGAGATAACTCACTGCGATATTTTCTCAGAATACTAGAGATAATAACC")


# === 2. Split into ~170bp windows ===
window_size = 180
windows = [seq[i:i+window_size] for i in range(0, len(seq), window_size)]
labels = [f"subunit_{i+1}" for i in range(len(windows))]

print("[INFO] Number of subunits:", len(windows))

# === 3. Define function to get percent identity ===
def percent_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    aln = alignments[0]
    matches = sum(a == b for a, b in zip(aln.seqA, aln.seqB) if a != '-' and b != '-')
    length = max(len(seq1), len(seq2))
    return 100 * matches / length

# === 4. Compute pairwise matrix ===
n = len(windows)
matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        pid = percent_identity(windows[i], windows[j])
        matrix[i, j] = pid

# === 5. Show matrix ===
df = pd.DataFrame(matrix, columns=labels, index=labels)
print(df.round(1))