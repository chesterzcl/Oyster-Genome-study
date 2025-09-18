import os
import re

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")

# === Configuration ===
EMAPPER_FILE = "eggnog_annotation.emapper.annotations"
GFF_FILE = "braker.gff3"
OUTPUT_GFF = "braker_with_name_and_description.gff3"
EVALUE_THRESHOLD = 1e-5

# === Step 1: Parse header from eggNOG file ===
print("Parsing eggNOG annotation...")
best_hits = {}

with open(EMAPPER_FILE) as f:
    # find header line
    for line in f:
        if line.startswith("query"):
            header = line.strip().split("\t")
            col_query = header.index("query")
            col_evalue = header.index("evalue")
            col_desc = header.index("Description")
            col_name = header.index("Preferred_name")
            break
    else:
        raise ValueError("No header line found in eggNOG file")

    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) <= max(col_desc, col_name):
            continue
        try:
            query = parts[col_query]
            evalue = float(parts[col_evalue])
            desc = parts[col_desc].strip()
            name = parts[col_name].strip()
        except (IndexError, ValueError):
            continue

        if (desc == "-" and name == "-") or evalue > EVALUE_THRESHOLD:
            continue

        if query not in best_hits or evalue < best_hits[query][0]:
            best_hits[query] = (evalue, name, desc)

print(f"Found {len(best_hits)} high-confidence annotations.")

# === Step 2: Canonicalize transcript IDs (e.g. g1.t1 -> g1) ===
def canonicalize(gene_id):
    return re.sub(r"\.t\d+$", "", gene_id)

gene_annot_map = {
    canonicalize(query): {"name": name, "desc": desc}
    for query, (_, name, desc) in best_hits.items()
}

# === Step 3: Annotate GFF ===
print("Annotating GFF...")
with open(GFF_FILE) as gff_in, open(OUTPUT_GFF, "w") as gff_out:
    for line in gff_in:
        if line.startswith("#") or "\t" not in line:
            gff_out.write(line)
            continue

        cols = line.strip().split("\t")
        attributes = cols[8]
        attr_dict = {}
        for attr in attributes.split(";"):
            if "=" in attr:
                k, v = attr.split("=", 1)
                attr_dict[k] = v

        gff_id = attr_dict.get("ID") or attr_dict.get("Parent")
        canonical_id = canonicalize(gff_id) if gff_id else None

        if canonical_id in gene_annot_map:
            name = gene_annot_map[canonical_id]["name"]
            desc = gene_annot_map[canonical_id]["desc"]
            if name and name != "-":
                attr_dict["Name"] = name
            if desc and desc != "-":
                attr_dict["Description"] = desc

        # Rebuild attributes
        cols[8] = ";".join(f"{k}={v}" for k, v in attr_dict.items())
        gff_out.write("\t".join(cols) + "\n")

print(f"Annotated GFF written to: {OUTPUT_GFF}")


BED_OUTPUT = "genes_function.bed"

print(f"Extracting gene features to BED: {BED_OUTPUT}")

with open(OUTPUT_GFF) as gff_in, open(BED_OUTPUT, "w") as bed_out:
    for line in gff_in:
        if line.startswith("#") or "\t" not in line:
            continue
        fields = line.strip().split("\t")
        feature_type = fields[2]
        if feature_type != "gene":
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1  # BED is 0-based
        end = fields[4]
        strand = fields[6]

        # Parse attributes to get ID and Name
        attributes = fields[8]
        attr_dict = {}
        for attr in attributes.split(";"):
            if "=" in attr:
                k, v = attr.split("=", 1)
                attr_dict[k] = v

        gene_id = attr_dict.get("ID", ".")
        gene_name = attr_dict.get("Name", "")

        label = gene_id

        bed_out.write(f"{chrom}\t{start}\t{end}\t{label}\t{gene_name}\t{strand}\n")

print(f"BED file written: {BED_OUTPUT}")