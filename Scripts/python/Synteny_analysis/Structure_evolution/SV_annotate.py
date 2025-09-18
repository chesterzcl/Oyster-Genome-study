import pandas as pd
import pybedtools
import os

# === Paths ===
dir = "/Users/lizicheng/Desktop/Data/dog_genome/oyster/sv/"
vcf_path = dir + "mantle_sv_filtered.vcf"
gff_path = dir + "annotated_genes.gff3"
sv_bed_path = dir + "sv_tmp.bed"

# === Step 1: Convert SV VCF to BED ===
def vcf_to_bed(vcf_file, output_bed):
    with open(vcf_file) as fin, open(output_bed, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1]) - 1  # BED is 0-based
            info = dict(item.split("=") for item in fields[7].split(";") if "=" in item)
            end = int(info.get("END", start + 1))  # fallback if no END
            svtype = info.get("SVTYPE", "NA")
            fout.write(f"{chrom}\t{start}\t{end}\t{svtype}\n")
    return output_bed

vcf_to_bed(vcf_path, sv_bed_path)

# === Step 2: Parse GFF3 and extract gene regions ===
gff_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
gff_df = pd.read_csv(gff_path, sep="\t", comment="#", names=gff_cols)

# Keep only features of interest
features_of_interest = ["gene", "exon", "intron", "CDS"]
gff_df = gff_df[gff_df["feature"].isin(features_of_interest)].copy()

# Extract gene ID or parent info
def extract_id(attr):
    for field in attr.split(";"):
        if field.startswith("ID=") or field.startswith("Parent="):
            return field.split("=")[1]
    return "NA"

gff_df["gene_id"] = gff_df["attributes"].apply(extract_id)
gff_df["start"] = gff_df["start"] - 1  # BED is 0-based

# Export GFF regions as BED
gff_bed_path = dir + "gff_regions.bed"
gff_df[["chrom", "start", "end", "feature", "gene_id"]].to_csv(gff_bed_path, sep="\t", index=False, header=False)

# === Step 3: Intersect SVs with GFF BED ===
sv_bed = pybedtools.BedTool(sv_bed_path)
gff_bed = pybedtools.BedTool(gff_bed_path)

overlap = sv_bed.intersect(gff_bed, wa=True, wb=True)

# === Step 4: Annotate ===
annotated = []
for line in overlap:
    f = str(line).strip().split("\t")
    sv_chr, sv_start, sv_end, sv_type = f[:4]
    feat_chr, feat_start, feat_end, region_type, gene_id = f[4:]

    annotated.append({
        "sv_chr": sv_chr,
        "sv_start": int(sv_start),
        "sv_end": int(sv_end),
        "sv_type": sv_type,
        "region_type": region_type,
        "gene_id": gene_id
    })

annotated_df = pd.DataFrame(annotated)

# Optional: group by SV and aggregate features
summary_df = (
    annotated_df
    .groupby(["sv_chr", "sv_start", "sv_end", "sv_type"])
    .agg(region_types=("region_type", lambda x: ",".join(sorted(set(x)))),
         genes=("gene_id", lambda x: ",".join(sorted(set(x)))))
    .reset_index()
)

# === Save results ===
annotated_df.to_csv(dir + "sv_detailed_annotation.tsv", sep="\t", index=False)
summary_df.to_csv(dir + "sv_annotation_summary.tsv", sep="\t", index=False)

print("Annotation complete. Files saved:")
print(" - sv_detailed_annotation.tsv")
print(" - sv_annotation_summary.tsv")

# Count unique SVs that overlap any gene-related feature
sv_overlapping_genes = annotated_df[annotated_df["region_type"] == "gene"]

n_overlap_gene = sv_overlapping_genes[["sv_chr", "sv_start", "sv_end", "sv_type"]].drop_duplicates().shape[0]
print(f"Number of SVs overlapping gene regions: {n_overlap_gene}")

# Count unique SVs per region type
region_counts = (
    annotated_df
    .groupby("region_type")[["sv_chr", "sv_start", "sv_end", "sv_type"]]
    .apply(lambda df: df.drop_duplicates().shape[0])
    .sort_values(ascending=False)
)

print("Number of SVs overlapping each region type:")
print(region_counts)