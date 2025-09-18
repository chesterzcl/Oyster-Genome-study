from collections import defaultdict

def parse_gff(gff_file):
    gene_transcripts = defaultdict(set)
    transcript_features = defaultdict(lambda: defaultdict(list))
    gene_coords = {}  # gene_id -> (chrom, start, end, strand)

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attributes = cols
            attr_dict = dict(item.split("=", 1) for item in attributes.split(";") if "=" in item)

            if feature_type == "gene":
                gene_id = attr_dict.get("ID")
                if gene_id:
                    gene_coords[gene_id] = (chrom, int(start)-1, int(end), strand)  # BED is 0-based

            elif feature_type == "mRNA":
                transcript_id = attr_dict.get("ID")
                parent_gene = attr_dict.get("Parent")
                if transcript_id and parent_gene:
                    gene_transcripts[parent_gene].add(transcript_id)

            elif feature_type in {"exon", "intron", "CDS"}:
                parent_id = attr_dict.get("Parent")
                if parent_id:
                    transcript_features[parent_id][feature_type].append((int(start), int(end)))

    return gene_transcripts, transcript_features, gene_coords

def detect_retrogenes(gene_transcripts, transcript_features):
    retrogenes = []

    for gene_id, transcripts in gene_transcripts.items():
        no_introns = all(
            not transcript_features[tid].get("intron")
            for tid in transcripts
        )
        if no_introns:
            retrogenes.append(gene_id)

    return retrogenes

def write_bed(retrogene_ids, gene_coords, output_bed):
    with open(output_bed, 'w') as out:
        for gene_id in retrogene_ids:
            if gene_id in gene_coords:
                chrom, start, end, strand = gene_coords[gene_id]
                out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")

# === Run ===
dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/sv/"
gff_path = dir+"annotated_genes.gff3"  # <-- replace this with your actual file path
output_bed = dir+"retrogenes.bed"     # <-- Output BED file


gene_transcripts, transcript_features, gene_coords = parse_gff(gff_path)
# === Get all gene IDs from gene_coords
all_genes = set(gene_coords.keys())
retrogenes = detect_retrogenes(gene_transcripts, transcript_features)

# === Determine non-retrogene gene IDs
retrogene_set = set(retrogenes)
non_retrogenes = list(all_genes - retrogene_set)
write_bed(retrogenes, gene_coords, output_bed)

# === Output file for non-retrogenes
nonretrogene_bed = dir + "non_retrogenes.bed"
write_bed(non_retrogenes, gene_coords, nonretrogene_bed)

print(f"Non-retrogene BED written to: {nonretrogene_bed}")
print(f"Non-retrogene count: {len(non_retrogenes)}")
print(f"Retrogene BED written to: {output_bed}")


print("Potential retrogenes (no introns in any isoform):")
for gid in retrogenes:
    print(gid)
