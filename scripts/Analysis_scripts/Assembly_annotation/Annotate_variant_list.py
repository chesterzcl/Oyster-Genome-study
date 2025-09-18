import sys
from collections import defaultdict


dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/"

# Inputs
vcf_file =dir+ "20cv_df2_hoh_biSNP_filtered_sig_dist.txt"
gff_file = dir+"annotated_genes.gff3"
output_file = dir+"snp_gene_overlap.tsv"

genes_by_chr = defaultdict(list)

with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        attr_field = fields[8]

        # Parse attributes
        attrs = {}
        for pair in attr_field.split(";"):
            if "=" in pair:
                key, value = pair.split("=", 1)
                attrs[key] = value

        gene_id = attrs.get("ID", ".")
        gene_name = attrs.get("Name", attrs.get("gene_name", "."))
        genes_by_chr[chrom].append((start, end, gene_id, gene_name))

# === PARSE VCF AND ANNOTATE ===
with open(vcf_file) as vcf, open(output_file, "w") as out:
    out.write("CHROM\tPOS\tREF\tALT\tGENE_ID\tGENE_NAME\n")
    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom, pos, ref, alts = fields[0], int(fields[1]), fields[3], fields[4]
        for alt in alts.split(","):  # Handle multi-allelics
            gene_id = "."
            gene_name = "."
            for start, end, gid, gname in genes_by_chr.get(chrom, []):
                if start <= pos <= end:
                    gene_id = gid
                    gene_name = gname
                    break
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene_id}\t{gene_name}\n")