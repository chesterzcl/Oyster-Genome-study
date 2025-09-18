import matplotlib.pyplot as plt
import os
import re
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from pycircos import Gcircle, Garc
from scipy.stats import mannwhitneyu
import collections
import os
import pandas as pd


spe1="mg"
# spe1="oe"
# spe1="ma"
spe2="cf"

os.chdir(f"/Users/lizicheng/Desktop/Data/dog_genome/oyster/{spe1.capitalize()}_CV_final")

col_file =f"{spe1.capitalize()}_vs_CV_final.collinearity"
gff_file =f"{spe1.capitalize()}_vs_CV_final.gff"
fai_file1="primary_dedup_chr_masked_hp_sealed.fa.fai"
fai_file2=f"{spe1.capitalize()}_clean.fa.fai"

#label1="C.virginica_3.0"
label1="xbMagGiga1.1"
# label1="ASM2561291v2"
# label1="xbOstEdul1.1"
label2="C.virginica_new"

# color1="#88b378"
color1="#F6A600"
color2="#B57EDC"

MIN_ANCHORS = 5

# your existing chromosome naming scheme
mg_chroms = [f"{spe1}{i}" for i in range(10, 0, -1)]
cf_chroms = [f"{spe2}{i}" for i in range(1, 11)]

label_map = {}

# Add M. gigas labels
for i in range(10, 0, -1):
    chrom_id = f"{spe1}{i}"
    label_map[chrom_id] = f"chr{i}"

# Add C. virginica labels
for i in range(1, 11):
    chrom_id = f"{spe2}{i}"
    label_map[chrom_id] = f"chr{i}"

# Print result to check
# print(label_map)


# Step 1: Parse GFF
def parse_fai_lengths(fai_path, prefix):
    """
    Returns {prefix+chrom_name: length}
    For example: prefix='cf' or 'mg'
    """
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            length = int(parts[1])
            # Add prefix to match naming scheme in your synteny
            lengths[prefix + chrom] = length
    return lengths

def parse_gff(gff_path):
    gene_pos = {}  # gene_name -> (scaffold, position)
    with open(gff_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue
            scaffold, gene, pos = parts
            pos = int(pos)
            gene_pos[gene] = (scaffold, pos)
    return gene_pos

def parse_collinearity_blocks(col_file):
    blocks = []
    block_id, orientation, chrA, chrB, block_pairs = None, None, None, None, []
    with open(col_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("## Alignment"):
                # If there was a previous block, save it
                if block_id is not None and block_pairs:
                    blocks.append((block_id, orientation, chrA, chrB, list(block_pairs)))
                # Parse the header
                m = re.match(r"## Alignment (\d+):.* ([^ ]+)&([^ ]+) (plus|minus)", line)
                if m:
                    block_id = int(m.group(1))
                    chrA, chrB, orientation = m.group(2), m.group(3), m.group(4)
                    block_pairs = []
                else:
                    # Malformed header, skip
                    block_id, orientation, chrA, chrB, block_pairs = None, None, None, None, []
            elif re.match(r"^\d+-\s*\d+:", line):
                # Example: 481-  0: g31251  gene-LOC105321635     5e-26
                parts = line.split()
                if len(parts) >= 4:
                    geneA = parts[2]
                    geneB = parts[3]
                    block_pairs.append((geneA, geneB))
        # Don't forget last block
        if block_id is not None and block_pairs:
            blocks.append((block_id, orientation, chrA, chrB, list(block_pairs)))
    return blocks

def block_ranges_from_pairs(blocks, gene_pos):
    block_ranges = []
    for block_id, orientation, chrA, chrB, pairs in blocks:
        posA = []
        posB = []
        for geneA, geneB in pairs:
            if geneA in gene_pos and geneB in gene_pos:
                cA, pA = gene_pos[geneA]
                cB, pB = gene_pos[geneB]
                if cA == chrA and cB == chrB:  # Strict match to block chrom
                    posA.append(pA)
                    posB.append(pB)
        if posA and posB:
            startA, endA = min(posA), max(posA)
            startB, endB = min(posB), max(posB)
            block_ranges.append((chrA, startA, endA, chrB, startB, endB, orientation))
    return block_ranges


# Step 3: Build chromosome layout (vertical bars)
def layout_chromosomes(chrom_lengths, x):
    y_offset = 0
    layout = {}
    spacing = 1000000  # 1Mb between bars
    for chrom, size in sorted(chrom_lengths.items()):
        layout[chrom] = (x, y_offset, y_offset + size)
        y_offset += size + spacing
    return layout


def plot_circos_blocks(gene_pos, chrom_lengths, block_ranges):
    from pycircos import Gcircle, Garc
    import matplotlib.pyplot as plt

    # --- Assign colors by genome ---
    def get_color(chrom):
        if chrom.startswith(spe1):
            return color1
        elif chrom.startswith(spe2):
            return color2
        else:
            return "#AAAAAA"

    # ---- Order chromosomes explicitly ----
    mg_chroms = [f"{spe1}{i}" for i in range(10, 0, -1)] 
    cv_chroms = [f"{spe2}{i}" for i in range(1, 11)]
    chrom_order = mg_chroms + cv_chroms

    circle = Gcircle()

    # Add arcs
    for chrom in chrom_order:
        if chrom in chrom_lengths:
            arc = Garc(
                arc_id=chrom,
                size=chrom_lengths[chrom],
                interspace=3,
                raxis_range=(850, 900),
                labelposition=60,
                label_visible=True,
                label=label_map.get(chrom, chrom),
                facecolor=get_color(chrom),
                edgecolor=get_color(chrom)
            )
            circle.add_garc(arc)

    circle.set_garcs()

    # Plot chords using precomputed block_ranges
    for chrA, startA, endA, chrB, startB, endB, orientation in block_ranges:
        color = "#4A90E2" if orientation == "plus" else "#D0021B"
        try:
            circle.chord_plot(
                (chrA, startA, endA, 800),
                (chrB, startB, endB, 800),
                facecolor=color,
                edgecolor="none"
            )
        except Exception as e:
            print(f"Skip block {chrA}:{startA}-{endA} <-> {chrB}:{startB}-{endB}: {e}")

    legend_handles = [
        Patch(color=color1, label=label1), 
        Patch(color=color2, label=label2),
        Patch(color="#4A90E2", label="Synteny: Same Strand"),
        Patch(color="#D0021B", label="Synteny: Inverted Strand"),
    ]
    plt.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(0.78, 0.93), fontsize=8,frameon=False)
    plt.savefig("Synteny_alignment_circos.png")


from collections import defaultdict

def merge_overlapping_blocks_with_orientation(block_ranges):
    """
    Input: block_ranges = list of (chrA, startA, endA, chrB, startB, endB, orientation)
    Output: merged_ranges = list of (chrA, merged_startA, merged_endA, chrB, merged_startB, merged_endB, orientation)
    """

    from collections import defaultdict

    # Group by (chrA, chrB, orientation)
    grouped = defaultdict(list)
    for br in block_ranges:
        key = (br[0], br[3], br[6])  # (chrA, chrB, orientation)
        grouped[key].append((br[1], br[2], br[4], br[5]))  # startA, endA, startB, endB

    merged_results = []

    for (chrA, chrB, orientation), intervals in grouped.items():
        # Sort by startA
        intervals.sort()

        current = list(intervals[0])

        for iv in intervals[1:]:
            # Check overlap in genome A coordinates
            if iv[0] <= current[1]:
                # Merge intervals in both A and B
                current[1] = max(current[1], iv[1])
                current[3] = max(current[3], iv[3])
                current[2] = min(current[2], iv[2])
            else:
                # Add completed merged interval
                merged_results.append((chrA, current[0], current[1], chrB, current[2], current[3], orientation))
                current = list(iv)

        # Add last interval
        merged_results.append((chrA, current[0], current[1], chrB, current[2], current[3], orientation))

    return merged_results


# Run
gene_pos = parse_gff(gff_file)
blocks=parse_collinearity_blocks(col_file)

# Extract block sizes
block_sizes = [len(b[4]) for b in blocks]

print(f"Total blocks: {len(blocks)}")
print(f"Block size range: {min(block_sizes)} - {max(block_sizes)}")

# Plot histogram
plt.figure(figsize=(8,5))
plt.hist(block_sizes, bins=50, color="steelblue", edgecolor="black")
plt.xlabel("Number of anchor pairs in block")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("block_size_distribution.png")
plt.close()

print("*** Wrote block_size_distribution.png")

cv_chrom_lengths = parse_fai_lengths(fai_file1, prefix="")
mg_chrom_lengths = parse_fai_lengths(fai_file2, prefix="")


####Filter by block size

filtered_blocks = [b for b in blocks if len(b[4]) >= MIN_ANCHORS]

print(f"Kept {len(filtered_blocks)} of {len(blocks)} total blocks.")

# Remap to match synteny naming
# mg_fai_map = {
#     f"chr{i}": f"mg{i}" for i in range(1, 11)
# }
mg_fai_map = {
    f"chr{i}": f"{spe1}{i}" for i in range(1, 11)
}
cf_fai_map = {
    f"HiC_scaffold_{i}": f"{spe2}{i}" for i in range(1, 11)
}

cv_chrom_lengths = {cf_fai_map.get(k, k): v for k, v in cv_chrom_lengths.items()}
mg_chrom_lengths = {mg_fai_map.get(k, k): v for k, v in mg_chrom_lengths.items()}

chrom_lengths = {}
chrom_lengths.update(cv_chrom_lengths)
chrom_lengths.update(mg_chrom_lengths)


block_ranges=block_ranges_from_pairs(filtered_blocks,gene_pos)
# Sum plus/minus alignment for each Mg chromosome
plus_lengths=defaultdict(int)
minus_lengths=defaultdict(int)

for chrA,startA,endA,chrB,startB,endB,orientation in block_ranges:
    length=abs(endA-startA)
    if orientation=='plus':
        plus_lengths[chrB]+=length
    elif orientation=='minus':
        minus_lengths[chrB]+=length

flip_mg_chroms=set()

all_mg_chroms_in_alignment=set(plus_lengths.keys())|set(minus_lengths.keys())

for chrom in all_mg_chroms_in_alignment:
    plus_total=plus_lengths.get(chrom,0)
    minus_total=minus_lengths.get(chrom,0)
    if minus_total>plus_total:
        flip_mg_chroms.add(chrom)
        print(f"Suggest flipping {chrom}: minus {minus_total} > plus {plus_total}")


# flip_mg_chroms.add('mg9')
# flip_mg_chroms.add('mg7')
# flip_mg_chroms.remove('mg1')
# flip_mg_chroms.remove('ma2')
# flip_mg_chroms.add('ma4')
# flip_mg_chroms.remove('ma5')
# flip_mg_chroms.remove('ma10')

print(flip_mg_chroms)
corrected_block_ranges = []

for chrA,startA,endA,chrB,startB,endB,orientation in block_ranges:
    if chrB in flip_mg_chroms:
        # Flip MG coordinates
        chrom_len=chrom_lengths[chrB]
        new_startB=chrom_len-endB
        new_endB=chrom_len-startB
        # Reverse orientation
        new_orientation='plus' if orientation=='minus' else 'minus'
    else:
        new_startB,new_endB=startB,endB
        new_orientation=orientation

    corrected_block_ranges.append(
        (chrA,startA,endA,chrB,new_startB,new_endB,new_orientation)
    )

merged_block_ranges=merge_overlapping_blocks_with_orientation(corrected_block_ranges)
merged_block_ranges=corrected_block_ranges

plot_circos_blocks(gene_pos,chrom_lengths,merged_block_ranges)


# Compute lengths in each genome
block_lengths_A = []
block_lengths_B = []
orientations = []
comparisons = []
chrA_list=[]

comparison_label = f"{spe1.upper()}_vs_CV"

for chrA, startA, endA, chrB, startB, endB, orientation in merged_block_ranges:
    lengthA = abs(endA-startA)
    lengthB = abs(endB-startB)
    block_lengths_A.append(lengthA)
    block_lengths_B.append(lengthB)
    orientations.append(orientation)
    comparisons.append(comparison_label)
    chrA_list.append(chrA)

print(f"Mean block length (CV): {np.mean(block_lengths_A):.1f} bp")
print(f"Mean block length ({spe1.capitalize()}): {np.mean(block_lengths_B):.1f} bp")

# Build DataFrame for export or plotting
df_synteny = pd.DataFrame({
    'block_length_cv': block_lengths_A,
    'block_length_spe1': block_lengths_B,
    'orientation': orientations,
    'comparison': comparisons,
    'cv_chromosome':chrA_list
})

# Export to TSV for later merge or plotting
df_synteny.to_csv(f"synteny_block_lengths_{spe1}_vs_cv.tsv", sep="\t", index=False)


total_span_A = sum(block_lengths_A)
total_span_B = sum(block_lengths_B)
print(f"Total syntenic span in {spe1.capitalize()}: {total_span_A/1e6:.2f} Mb")
print(f"Total syntenic span in CV: {total_span_B/1e6:.2f} Mb")


total_genome_size_A = sum([v for k,v in chrom_lengths.items() if k.startswith(spe1)])
total_genome_size_B = sum([v for k,v in chrom_lengths.items() if k.startswith(spe2)])

print(f"Proportion syntenic in {spe1.capitalize()}: {total_span_A/total_genome_size_A:.2%}")
print(f"Proportion syntenic in CV: {total_span_B/total_genome_size_B:.2%}")

from collections import Counter
orientation_counts = Counter(orientations)
print("Orientation counts:", orientation_counts)

plt.figure(figsize=(8,5))
plt.hist(np.log10(block_lengths_A), bins=30, alpha=0.7, label=label1)
plt.hist(np.log10(block_lengths_B), bins=30, alpha=0.7, label=label2)
plt.xlabel('log10(Block Length) [bp]')
plt.ylabel('Number of Blocks')
plt.legend()
plt.savefig('Synteny_block_length.png')
# plt.show()

df_blocks = pd.DataFrame(merged_block_ranges, columns=['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'orientation'])
df_blocks['lengthA'] = df_blocks['endA'] - df_blocks['startA']
df_blocks['lengthB'] = df_blocks['endB'] - df_blocks['startB']

df_blocks.to_csv('synteny_block_summary.csv', index=False)
print("Saved synteny_block_summary.csv")


# Split corrected_block_ranges by C. virginica chromosome
blocks_by_cf = defaultdict(list)
for block in merged_block_ranges:
    chrA, startA, endA, chrB, startB, endB, orientation = block
    blocks_by_cf[chrA].append(block)

# Sort each Cv chromosome's blocks by Cv coordinates
for cf_chrom in blocks_by_cf:
    blocks_by_cf[cf_chrom].sort(key=lambda x: x[1])


orientation_flips = []

for cf_chrom in blocks_by_cf:
    blocks = blocks_by_cf[cf_chrom]
    for i in range(1, len(blocks)):
        prev = blocks[i-1]
        curr = blocks[i]
        _, _, _, prev_mg, _, _, prev_orient = prev
        _, _, _, curr_mg, _, _, curr_orient = curr

        # Only consider same MG chromosome
        if prev_mg == curr_mg and prev_orient != curr_orient:
            orientation_flips.append((cf_chrom, prev_mg, prev_orient, curr_orient))

# Summarize flips
from collections import Counter
flip_counts = Counter((cf, mg) for cf, mg, _, _ in orientation_flips)

# print("\nOrientation flips per chromosome pair:")
# for (cf, mg), count in flip_counts.items():
#     print(f"{cf} vs {mg}: {count} flips")

from scipy.stats import spearmanr


spearman_stats = {}   # (cf, mg) -> (rho, pvalue)

for cf_chrom in blocks_by_cf:
    blocks = blocks_by_cf[cf_chrom]
    mg_chrom_set = set(b[3] for b in blocks)

    for mg_chrom in mg_chrom_set:
        pair_blocks = [b for b in blocks if b[3] == mg_chrom]
        if len(pair_blocks) < 5:
            continue

        cf_starts = [b[1] for b in pair_blocks]
        mg_starts = [b[4] for b in pair_blocks]

        rho, p = spearmanr(cf_starts, mg_starts)
        spearman_stats[(cf_chrom, mg_chrom)] = (rho, p)

        # print(f"{cf_chrom} vs {mg_chrom}: Spearman rho={rho:.2f} (p={p:.3g})")


# Store summed inverted lengths per cf-mg pair
inverted_lengths=defaultdict(int)

for b in corrected_block_ranges:
    cf_chrom,mg_chrom=b[0],b[3]
    orientation=b[6]
    block_length=abs(b[2]-b[1])
    
    if orientation=='minus':
        inverted_lengths[(cf_chrom,mg_chrom)]+=block_length

# print("\nTotal inverted length per chromosome pair:")
# for (cf, mg), length in inverted_lengths.items():
#     print(f"{cf} vs {mg}: {length/1e6:.2f} Mb")

# Count total and minus-strand blocks
num_blocks=defaultdict(int)
num_flipped_blocks=defaultdict(int)

for b in corrected_block_ranges:
    cf_chr,mg_chr=b[0],b[3]
    orientation=b[6]
    
    num_blocks[(cf_chr,mg_chr)]+=1
    if orientation=='minus':
        num_flipped_blocks[(cf_chr, mg_chr)]+=1


# Make sure we get all pairs seen in any of the stats
all_pairs=(
    set(flip_counts.keys())|
    set(inverted_lengths.keys())|
    set(spearman_stats.keys())|
    set(num_blocks.keys())
)


def count_breakpoints(order):
    count=0
    for i in range(len(order)-1):
        if order[i+1]!=order[i]+1:
            count+=1
    return count

import csv

all_pairs=(
    set(flip_counts.keys())|
    set(inverted_lengths.keys())|
    set(spearman_stats.keys())|
    set(num_blocks.keys())
)


with open("global_rearrangement_summary.csv", "w", newline="") as csvfile:
    writer=csv.writer(csvfile)
    writer.writerow([
        f"CV chromosome", f"{spe1.capitalize()} chromosome", 
        "# Syntenic Blocks", "Reversed Block Length (Mb)", 
        "% Blocks Reversed","# Orientation Flips",
        "Spearman ρ", "p-value",
        # "Breakpoints", "#Breakpoints",
    ])
    
    for (cf, mg) in sorted(all_pairs):
        total_blocks=num_blocks.get((cf, mg), 0)
        flipped_blocks=num_flipped_blocks.get((cf, mg), 0)
        percent_flipped=100.0*flipped_blocks/total_blocks if total_blocks else 0.0
        flips=flip_counts.get((cf,mg),0)
        flip_len_Mb=inverted_lengths.get((cf,mg),0)/1e6

        # Spearman stats
        rho,pval=spearman_stats.get((cf,mg),(None,None))
        rho_str=f"{rho:.3f}" if rho is not None else ""
        pval_str=f"{pval:.3g}" if pval is not None else ""

        # ------------------------------
        # Compute breakpoints
        # ------------------------------
        # Get the blocks for this pair
        pair_blocks = [
            b for b in blocks_by_cf.get(cf, [])
            if b[3] == mg
        ]

        # Require at least 5 blocks for this metric
        if len(pair_blocks) >= 5:
            # They are already sorted by cf coordinates
            mg_positions = [b[4] for b in pair_blocks]
            # Get MG rank order
            mg_ranks = sorted(range(len(mg_positions)), key=lambda i: mg_positions[i])
            mg_order = [mg_ranks.index(i)+1 for i in range(len(mg_positions))]
            # Count breakpoints
            breakpoints = count_breakpoints(mg_order)
            max_breakpoints = len(mg_order) - 1
            breakpoints_percent = 100.0 * breakpoints / max_breakpoints if max_breakpoints > 0 else 0.0
        else:
            breakpoints = ""
            breakpoints_percent = ""

        # ------------------------------
        # Write row
        # ------------------------------
        writer.writerow([
            cf, mg,
            total_blocks,f"{flip_len_Mb:.2f}",
            f"{percent_flipped:.1f}",flips,
            rho_str, pval_str,
            # breakpoints, f"{breakpoints_percent:.1f}" if breakpoints_percent != "" else ""
        ])




# blocks = parse_collinearity_blocks(col_file)

block_id_to_pairs = {
    b[0]: b[4]
    for b in filtered_blocks
}

# -------------------------------------------------------------
# 2. Compute block ranges with *original* orientation
# -------------------------------------------------------------
block_ranges = []
for block_id, orientation, chrA, chrB, pairs in filtered_blocks:
    posA = []
    posB = []
    for geneA, geneB in pairs:
        if geneA in gene_pos and geneB in gene_pos:
            cA, pA = gene_pos[geneA]
            cB, pB = gene_pos[geneB]
            if cA == chrA and cB == chrB:
                posA.append(pA)
                posB.append(pB)
    if posA and posB:
        startA, endA = min(posA), max(posA)
        startB, endB = min(posB), max(posB)
        block_ranges.append((block_id, chrA, startA, endA, chrB, startB, endB, orientation))


# -------------------------------------------------------------
# 3. Apply orientation correction (flip MG coords if needed)
# -------------------------------------------------------------
corrected_block_ranges = []
for block_id, chrA, startA, endA, chrB, startB, endB, orientation in block_ranges:
    if chrB in flip_mg_chroms:
        chrom_len = chrom_lengths[chrB]
        new_startB = chrom_len - endB
        new_endB = chrom_len - startB
        new_orientation = 'plus' if orientation == 'minus' else 'minus'
    else:
        new_startB, new_endB = startB, endB
        new_orientation = orientation

    corrected_block_ranges.append((
        block_id, chrA, startA, endA, chrB, new_startB, new_endB, new_orientation
    ))


# -------------------------------------------------------------
# 4. Extract ANCHOR gene pairs (with orientation correction)
# -------------------------------------------------------------
gene_pair_rows = []
for block_id, chrA, startA, endA, chrB, startB, endB, orientation in corrected_block_ranges:
    anchor_pairs = block_id_to_pairs.get(block_id, [])
    for geneA, geneB in anchor_pairs:
        if geneA in gene_pos and geneB in gene_pos:
            scaffoldA, posA = gene_pos[geneA]
            scaffoldB, posB = gene_pos[geneB]

            if scaffoldA != chrA or scaffoldB != chrB:
                continue

            # Flip MG gene position if needed
            if chrB in flip_mg_chroms:
                chrom_len = chrom_lengths[chrB]
                posB = chrom_len - posB

            gene_pair_rows.append((
                block_id, orientation,
                chrA, posA, geneA,
                chrB, posB, geneB
            ))

# -------------------------------------------------------------
# 5. Write ANCHOR genes to DataFrame
# -------------------------------------------------------------
import pandas as pd
df_pairs = pd.DataFrame(
    gene_pair_rows,
    columns=["block_id", "orientation", "chrA", "posA", "geneA", "chrB", "posB", "geneB"]
)
print(f"Extracted {len(df_pairs)} anchor gene pairs")
df_pairs.to_csv("synteny_gene_pairs.csv", index=False)

# -------------------------------------------------------------
# 6. Determine NON-ANCHOR genes in block ranges
# -------------------------------------------------------------
anchor_set_A = set(zip(df_pairs["chrA"], df_pairs["geneA"]))
anchor_set_B = set(zip(df_pairs["chrB"], df_pairs["geneB"]))

non_anchor_records = []
for block_id, chrA, startA, endA, chrB, startB, endB, orientation in corrected_block_ranges:
    # Genes on Cv
    for gene, pos in [(g, p) for g, (sc, p) in gene_pos.items() if sc == chrA and startA <= p <= endA]:
        if (chrA, gene) not in anchor_set_A:
            non_anchor_records.append([
                block_id, chrA, gene, pos, "", "", "", orientation, "non-anchor"
            ])

    # Genes on Mg
    for gene, pos in [(g, p) for g, (sc, p) in gene_pos.items() if sc == chrB and startB <= p <= endB]:
        if (chrB, gene) not in anchor_set_B:
            non_anchor_records.append([
                block_id, "", "", "", chrB, gene, pos, orientation, "non-anchor"
            ])

print(f"Found {len(non_anchor_records)} non-anchor genes")

# -------------------------------------------------------------
# 7. Combine and write ANNOTATED table
# -------------------------------------------------------------
combined_records = []
for _, row in df_pairs.iterrows():
    combined_records.append([
        row["block_id"],
        row["chrA"], row["geneA"], row["posA"],
        row["chrB"], row["geneB"], row["posB"],
        row["orientation"],
        "anchor"
    ])
combined_records.extend(non_anchor_records)

import csv

with open("synteny_genes_annotated.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "block_id",
        "chrA", "geneA", "posA",
        "chrB", "geneB", "posB",
        "orientation",
        "type"
    ])
    writer.writerows(combined_records)

print("Wrote synteny_genes_annotated.csv with anchors and non-anchors")


# =====================================
# 0. Assumes you already loaded combined_records
# =====================================
print(f"Loaded {len(combined_records)} records from synteny_genes_annotated.csv")

# -------------------------------------
# Prepare gene lookup tables for CF and MG
# -------------------------------------
gene_list_cf = defaultdict(list)
gene_list_mg = defaultdict(list)

for row in combined_records:
    _, chrA, geneA, posA, chrB, geneB, posB, orientation, _type = row
    if geneA:
        gene_list_cf[chrA].append((int(posA), geneA))
    if geneB:
        gene_list_mg[chrB].append((int(posB), geneB))

for chrom in gene_list_cf:
    gene_list_cf[chrom].sort()
for chrom in gene_list_mg:
    gene_list_mg[chrom].sort()

print("*** Prepared gene lookup tables.")
# print(gene_list_cf)
# =====================================
# 1. Filter only anchor records
# =====================================
anchor_records = [row for row in combined_records if row[-1] == "anchor"]
print(f"Found {len(anchor_records)} anchor records")

# Group by cf chromosome
anchors_by_cf = defaultdict(list)
for row in anchor_records:
    block_id, chrA, geneA, posA, chrB, geneB, posB, orientation, _type = row
    anchors_by_cf[chrA].append({
        "block_id": int(block_id),
        "geneA": geneA,
        "posA": int(posA),
        "chrB": chrB,
        "geneB": geneB,
        "posB": int(posB),
        "orientation": orientation
    })

# =====================================
# 2. Compute raw adjacent gap distributions
# =====================================
cf_gaps_all = []
mg_gaps_all = []
gap_records_basic = []

for cf_chr, anchors in anchors_by_cf.items():
    if len(anchors) < 2:
        continue
    anchors.sort(key=lambda x: x["posA"])

    # Remove exact duplicates
    deduped = []
    last_key = None
    for a in anchors:
        key = (a["geneA"], a["posA"])
        if key != last_key:
            deduped.append(a)
        last_key = key
    anchors = deduped

    for i in range(len(anchors) - 1):
        a1 = anchors[i]
        a2 = anchors[i+1]

        # Skip "self" artifacts
        if (a1["geneA"] == a2["geneA"]) and (a1["posA"] == a2["posA"]):
            continue

        cf_gap = a2["posA"] - a1["posA"]
        mg_gap = abs(a2["posB"] - a1["posB"])
        orientation_change = a1["orientation"] != a2["orientation"]
        boundary_label = "within_block" if a1["block_id"] == a2["block_id"] else "boundary"

        cf_gaps_all.append(cf_gap)
        mg_gaps_all.append(mg_gap)
        gap_records_basic.append([
            cf_chr,
            a1["block_id"], a2["block_id"],
            a1["geneA"], a1["posA"],
            a2["geneA"], a2["posA"],
            cf_gap, mg_gap,
            a1["chrB"],
            a1["geneB"], a1["posB"],
            a2["geneB"], a2["posB"],
            orientation_change,
            boundary_label
        ])

print(f"Collected {len(cf_gaps_all)} adjacent anchor pairs.")

# =====================================
# 3. Plot gap size histograms
# =====================================
fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(np.log10([g for g in cf_gaps_all if g > 0]), bins=100, alpha=1, label=label2,color=color2)
ax.hist(np.log10([g for g in mg_gaps_all if g > 0]), bins=100, alpha=0.4, label=label1,color=color1)
# sns.kdeplot(np.log10([g for g in cf_gaps_all if g > 0]), label=label2, color=color2, fill=False, alpha=0.5)
# sns.kdeplot(np.log10([g for g in mg_gaps_all if g > 0]), label=label1, color=color1, fill=False, alpha=0.5)
# plt.xlabel("log10(gap size) [bp]",fontsize=18)
plt.ylabel("Count",fontsize=18)
ax.legend(fontsize=18,frameon=False)
ax.set_xlim([2.8,8])
# ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(labelsize=18)
ax.set_xticks([4,6,8])       # Remove tick marks
ax.set_xticklabels([])  # Remove tick labels
plt.tight_layout()
plt.savefig("gap_size_distribution_combined.png")
plt.close()
print("Saved gap_size_distribution_combined.png")

# =====================================
# 4. Compute robust thresholds (99th percentile)
# =====================================
within_cf_gaps = [r[7] for r in gap_records_basic if r[-1] == "within_block" and not r[-2]]
within_mg_gaps = [r[8] for r in gap_records_basic if r[-1] == "within_block" and not r[-2]]

threshold_cf = np.percentile(within_cf_gaps, 99)
threshold_mg = np.percentile(within_mg_gaps, 99)

print(f"Thresholds chosen at 99th percentile of within_block distribution:")
print(f"  cf_gap threshold: {threshold_cf:.0f} bp")
print(f"  {spe1}_gap threshold: {threshold_mg:.0f} bp")

# Plot thresholds on within-block histograms
plt.figure(figsize=(6, 5))
plt.hist(np.log10([g for g in within_cf_gaps if g > 0]), bins=50, alpha=0.5, label=label2,color=color2)
plt.hist(np.log10([g for g in within_mg_gaps if g > 0]), bins=50, alpha=0.5, label=label1,color=color1)
plt.axvline(np.log10(threshold_cf), color=color2, linestyle='--')
plt.axvline(np.log10(threshold_mg), color=color1, linestyle='--')
plt.legend()
plt.xlabel('log10(gap size)')
plt.ylabel('Count')
plt.xlim(2,7)
plt.savefig('gap_distribution_99pct_threshold.png')
plt.close()
print("Saved gap_distribution_99pct_threshold.png")

# =====================================
# 5. Build anchor lookup for counting
# =====================================
anchors_by_mg = defaultdict(list)
for row in anchor_records:
    anchors_by_mg[row[4]].append({
        "geneB": row[5],
        "posB": int(row[6])
    })

for chrom in anchors_by_cf:
    anchors_by_cf[chrom].sort(key=lambda x: x["posA"])
for chrom in anchors_by_mg:
    anchors_by_mg[chrom].sort(key=lambda x: x["posB"])

print(f"*** Prepared sorted anchor lists for both cf and {spe1}.")

# =====================================
# 6. Classify and count anchors
# =====================================
def has_foreign_anchor_between(anchors_df, chrB, posB1, posB2, block_id1, block_id2):
    """
    Check if any anchor gene between posB1 and posB2 has a block_id 
    that is not block_id1 or block_id2 on chrB.

    Parameters:
        anchors_df (pd.DataFrame): Must contain ['block_id', 'chrB', 'posB']
        chrB (str): Chromosome name
        posB1 (int): Position of first anchor gene
        posB2 (int): Position of second anchor gene
        block_id1 (int or str): block_id of the first flanking anchor gene
        block_id2 (int or str): block_id of the second flanking anchor gene

    Returns:
        bool: True if foreign anchor(s) exist, else False
        pd.DataFrame: Details of foreign anchors
    """
    start, end = sorted([posB1, posB2])

    # Subset to anchor genes in range and chromosome
    region_genes = anchors_df[
        (anchors_df["chrB"] == chrB) &
        (anchors_df["posB"] > start) &
        (anchors_df["posB"] < end)
    ]

    # Keep genes with block_id not matching either block_id1 or block_id2
    foreign_genes = region_genes[
        ~region_genes["block_id"].isin([block_id1, block_id2])
    ]

    return not foreign_genes.empty, foreign_genes

gap_records = []

for cf_chr, anchors in anchors_by_cf.items():
    if len(anchors) < 2:
        continue
    anchors.sort(key=lambda x: x["posA"])

    for i in range(len(anchors) - 1):
        a1 = anchors[i]
        a2 = anchors[i+1]

        if (a1["geneA"] == a2["geneA"]) and (a1["posA"] == a2["posA"]):
            continue

        b1=a1['block_id']
        b2=a2['block_id']
        cf_gap = a2["posA"] - a1["posA"]
        mg_gap = abs(a2["posB"] - a1["posB"])
        orientation_change = a1["orientation"] != a2["orientation"]
        boundary_label = "within_block" if a1["block_id"] == a2["block_id"] else "boundary"

        # Event classification
        # First, detect complex cases involving both genomes
        if cf_gap > threshold_cf and mg_gap > threshold_mg:
            event_type = "complex_sv"

        # Now check clean inversion
        elif orientation_change:
            if cf_gap <= threshold_cf and mg_gap <= threshold_mg:
                event_type = "inversion"
            else:
                event_type = "complex_sv"

        # Insertion biases
        elif cf_gap > threshold_cf:
            event_type = "insertion_cv"

        elif mg_gap > threshold_mg:
            found, foreign = has_foreign_anchor_between(
                df_pairs,
                chrB=a1["chrB"],
                posB1=a1["posB"],
                posB2=a2["posB"],
                block_id1=b1,
                block_id2=b2  # Same block, e.g., consecutive anchors
            )
            if found:
                # print(foreign)
                event_type=f"insertion_{spe1}_foreign"
            else:
                event_type = f"insertion_{spe1}"

        # Default
        else:
            event_type = "trivial"

        # Count anchors (inclusive)
        cf_anchor_count = sum(
            a1["posA"] <= x["posA"] <= a2["posA"]
            for x in anchors_by_cf.get(cf_chr, [])
        )
        mg_anchor_count = sum(
            min(a1["posB"], a2["posB"]) <= x["posB"] <= max(a1["posB"], a2["posB"])
            for x in anchors_by_mg.get(a1["chrB"], [])
        )

        gap_records.append([
            cf_chr,
            a1["block_id"], a2["block_id"],
            a1["geneA"], a1["posA"],
            a2["geneA"], a2["posA"],
            cf_gap, mg_gap,
            a1["chrB"],
            a1["geneB"], a1["posB"],
            a2["geneB"], a2["posB"],
            orientation_change,
            boundary_label,
            event_type,
            cf_anchor_count,
            mg_anchor_count
        ])

print(f"Identified {len(gap_records)} adjacent anchor pairs with classification and anchor counts.")

# =====================================
# 7. Save to CSV
# =====================================
with open("gap_records_analysis.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "cf_chr", "block_id1", "block_id2",
        "cf_gene1", "cf_pos1", "cf_gene2", "cf_pos2",
        "cf_gap", f"{spe1}_gap",
        f"{spe1}_chr", f"{spe1}_gene1", f"{spe1}_pos1", f"{spe1}_gene2", f"{spe1}_pos2",
        "orientation_change",
        "boundary_label",
        "event_type",
        "cf_anchor_count",
        f"{spe1}_anchor_count"
    ])
    writer.writerows(gap_records)

print("*** Saved gap_records_analysis.csv")

# =====================================
# 8. Analyze and Plot
# =====================================
df_gaps = pd.DataFrame(gap_records, columns=[
    "cf_chr", "block_id1", "block_id2",
    "cf_gene1", "cf_pos1", "cf_gene2", "cf_pos2",
    "cf_gap", f"{spe1}_gap",
    f"{spe1}_chr", f"{spe1}_gene1", f"{spe1}_pos1", f"{spe1}_gene2", f"{spe1}_pos2",
    "orientation_change",
    "boundary_label",
    "event_type",
    "cf_anchor_count",
    f"{spe1}_anchor_count"
])

# Reshape to long format
long_gaps = pd.melt(
    df_gaps,
    id_vars=['event_type'],
    value_vars=['cf_gap', f"{spe1}_gap"],
    var_name='gap_type',
    value_name='gap_value'
)

# Remove zero or negative gaps if any
long_gaps = long_gaps[long_gaps['gap_value'] > 0]

long_gaps = long_gaps[long_gaps["event_type"] != "trivial"]
custom_order = [
    "inversion",
    "insertion_cv",
    f"insertion_{spe1}",
    f"insertion_{spe1}_foreign",
    "complex_sv",
]

# Log transform if desired
long_gaps['log_gap'] = np.log10(long_gaps['gap_value'])

custom_palette = {
    "cf_gap": color2,   
    f"{spe1}_gap": color1    
}

#====test
# Filter data
sub_df = long_gaps[(long_gaps["gap_type"] == f"{spe1}_gap") & 
                   (long_gaps["event_type"].isin([f"insertion_{spe1}", f"insertion_{spe1}_foreign"]))]

# Separate groups
group_mg = sub_df[sub_df["event_type"] == f"insertion_{spe1}"]["gap_value"]
group_foreign = sub_df[sub_df["event_type"] == f"insertion_{spe1}_foreign"]["gap_value"]

# Mann–Whitney U Test (two-sided by default)
stat, p = mannwhitneyu( group_foreign,group_mg, alternative="greater")
print(f"Mann–Whitney U test statistic: {stat:.3f}, p-value: {p:.4e}")

# ============================BoxPlot============================
plt.figure(figsize=(6, 10))
ax = sns.violinplot(
    y="event_type",
    x="log_gap",
    hue="gap_type",
    data=long_gaps,
    order=custom_order,
    # showfliers=False,
    width=0.8,
    palette=custom_palette
)
ax.set_ylabel("Event Type", fontsize=16)
ax.set_xlabel("log10(Gap Size) [bp]", fontsize=16)
ax.set_xlim([2.8,8])
# ============================BoxPlot============================

# Grab original colored legend handles
handles, labels = ax.get_legend_handles_labels()
# Replace only the text labels while keeping colors
ax.legend(
    handles=handles,
    labels=[label2,label1],
    # title="Gap Type",
    loc='upper left',
    frameon=False,
    fontsize=14,
    title_fontsize=14,
    bbox_to_anchor=(0, 1.15)
)
ax.set_yticklabels([
    "Inversion",
    "Deletion",
    "Local Insertion",
    "Foreign Insertion",
    "Complex SV"
], fontsize=18)
ax.tick_params(labelsize=18)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('none')  # No ticks on y-axis
ax.xaxis.set_ticks_position('bottom')  # Only bottom ticks

plt.tight_layout()
plt.savefig("grouped_gap_by_event_type_box.png")
plt.close()

print("*** Saved grouped boxplot as grouped_gap_boxplot_by_event_type.png ***")

# ============================BarPlot============================
count_df = (
    long_gaps
    .groupby(['event_type'])
    .size()
    .reset_index(name='event_count')
)
print(count_df)
plt.figure(figsize=(3, 10))
ax = sns.barplot(
    data=count_df,
    y="event_type",
    x="event_count",
    # hue="gap_type",
    order=custom_order,
    width=0.3,
    alpha=0.7
    # palette=custom_palette
)

# ax.set_yticklabels([
#     "Inversion",
#     "Deletion",
#     "Local Insertion",
#     "Foreign Insertion",
#     "Complex SV"
# ], fontsize=14)
ax.tick_params(labelsize=18)
ax.tick_params(axis='y', left=False, labelleft=False)
ax.yaxis.set_ticks_position('none') 
ax.set_ylabel("")
ax.set_xlabel("Event Count",fontsize=18)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig("grouped_gap_by_event_type_bar.png")
plt.close()
# ============================BarPlot============================

# =====================================
# Mann-Whitney U test
# =====================================
within_cf = df_gaps[df_gaps["boundary_label"]=="within_block"]["cf_gap"]
boundary_cf = df_gaps[df_gaps["boundary_label"]=="boundary"]["cf_gap"]
u_cf, p_cf = mannwhitneyu(within_cf, boundary_cf, alternative="two-sided")

within_mg = df_gaps[df_gaps["boundary_label"]=="within_block"][f"{spe1}_gap"]
boundary_mg = df_gaps[df_gaps["boundary_label"]=="boundary"][f"{spe1}_gap"]
u_mg, p_mg = mannwhitneyu(within_mg, boundary_mg, alternative="two-sided")

print("\n=== Statistical Test Results ===")
print(f"cf_gap Mann-Whitney U: U={u_cf}, p-value={p_cf}")
print(f"{spe1}_gap Mann-Whitney U: U={u_mg}, p-value={p_mg}")

print("\n*** DONE ***")


# =====================================
# Extracting gap genes with event-type-specific rules
# =====================================
# CONFIG
WINDOW_BP = 100_000
BUFFER_BP = 10_000


flanking_records = []

print("*** Extracting genes for boundary rows with non-trivial event types ***")

boundary_rows = [
    row for row in gap_records
    if row[15] == "boundary" and row[16] != "trivial"
]

print(f"Found {len(boundary_rows)} boundary rows (non-trivial)")

for row in boundary_rows:
    cf_chr = row[0]
    cf_pos1 = int(row[4])
    cf_pos2 = int(row[6])
    mg_chr = row[9]
    mg_pos1 = int(row[11])
    mg_pos2 = int(row[13])
    orientation_change = row[14]
    boundary_label = row[15]
    event_type = row[16]
    cf_gap = int(row[7])
    mg_gap = int(row[8])
    cf_genes = gene_list_cf.get(cf_chr, [])
    mg_genes = gene_list_mg.get(mg_chr, [])

    # Initialize
    cf_flank = []
    mg_flank = []

    if event_type in ["insertion_cv"]:
        # Inside the gap region in CF genome
        cf_region_start = min(cf_pos1, cf_pos2) - BUFFER_BP
        cf_region_end = max(cf_pos1, cf_pos2) + BUFFER_BP
        cf_flank = [(pos, gene) for pos, gene in cf_genes if cf_region_start <= pos <= cf_region_end]

    elif event_type in [f"insertion_{spe1}"]:
        # Inside the gap region in MG genome
        mg_region_start = min(mg_pos1, mg_pos2) - BUFFER_BP
        mg_region_end = max(mg_pos1, mg_pos2) + BUFFER_BP
        mg_flank = [(pos, gene) for pos, gene in mg_genes if mg_region_start <= pos <= mg_region_end]

    elif event_type in ["divergent_break", "inversion", "complex_sv"]:
        # 100kb windows around BOTH breakends in BOTH genomes
        cf_regions = [
            (cf_pos1 - WINDOW_BP, cf_pos1 + WINDOW_BP),
            (cf_pos2 - WINDOW_BP, cf_pos2 + WINDOW_BP)
        ]
        mg_regions = [
            (mg_pos1 - WINDOW_BP, mg_pos1 + WINDOW_BP),
            (mg_pos2 - WINDOW_BP, mg_pos2 + WINDOW_BP)
        ]

        for start, end in cf_regions:
            cf_flank.extend([(pos, gene) for pos, gene in cf_genes if start <= pos <= end])

        for start, end in mg_regions:
            mg_flank.extend([(pos, gene) for pos, gene in mg_genes if start <= pos <= end])

    else:
        # Other types are skipped
        continue

    # Skip if both sides empty
    if not cf_flank and not mg_flank:
        continue

    flanking_records.append([
        cf_chr, cf_pos1, cf_pos2,
        mg_chr, mg_pos1, mg_pos2,
        orientation_change, boundary_label, event_type,
        cf_gap, mg_gap,
        cf_flank,
        mg_flank
    ])


with open("breakpoint_flanking_genes.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        "cf_chr", "cf_pos1", "cf_pos2",
        f"{spe1}_chr", f"{spe1}_pos1", f"{spe1}_pos2",
        "orientation_change", "boundary_label", "event_type",
        "cf_gap", f"{spe1}_gap",
        "cf_flanking_genes", f"{spe1}_flanking_genes"
    ])
    
    for rec in flanking_records:
        (
            cf_chr, cf_pos1, cf_pos2,
            mg_chr, mg_pos1, mg_pos2,
            orientation_change, boundary_label, event_type,
            cf_gap, mg_gap,
            cf_flank, mg_flank
        ) = rec

        cf_flank_str = ";".join([f"{gene}@{pos}" for pos, gene in cf_flank])
        mg_flank_str = ";".join([f"{gene}@{pos}" for pos, gene in mg_flank])

        writer.writerow([
            cf_chr, cf_pos1, cf_pos2,
            mg_chr, mg_pos1, mg_pos2,
            orientation_change, boundary_label, event_type,
            cf_gap, mg_gap,
            cf_flank_str, mg_flank_str
        ])

print("*** DONE *** Wrote breakpoint_flanking_genes.csv")



