dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/"

import re

# Configuration
EMAPPER_FILE = "eggnog_annotation.emapper.annotations"
GFF_FILE = "braker.gff3"
OUTPUT_GFF = "braker_with_symbols.gff3"
EVALUE_THRESHOLD = 1e-5
EMAPPER_QUERY_COL = 0
EMAPPER_EVALUE_COL = 2
EMAPPER_SYMBOL_COL = 8  # "Preferred_name"

# Step 1: Parse emapper annotations (keep best hit per query with symbol and good E-value)
print("Parsing emapper file...")
best_hits = {}

with open(dir+EMAPPER_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) <= EMAPPER_SYMBOL_COL:
            continue

        query = parts[EMAPPER_QUERY_COL]
        symbol = parts[EMAPPER_SYMBOL_COL]
        try:
            evalue = float(parts[EMAPPER_EVALUE_COL])
        except ValueError:
            continue

        if symbol == "-" or evalue > EVALUE_THRESHOLD:
            continue

        if query not in best_hits or evalue < best_hits[query][0]:
            best_hits[query] = (evalue, symbol)

print(f"Found {len(best_hits)} high-confidence annotated transcripts.")

# Step 2: Build canonical map (strip .t1, .t2 suffix)
def canonicalize(gene_id):
    return re.sub(r"\.t\d+$", "", gene_id)

gene_symbol_map = {
    canonicalize(query): symbol for query, (evalue, symbol) in best_hits.items()
}

# Step 3: Annotate the GFF3 file
print("Annotating GFF file...")
with open(dir+GFF_FILE) as gff_in, open(dir+OUTPUT_GFF, "w") as gff_out:
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

        if canonical_id in gene_symbol_map:
            attr_dict["Name"] = gene_symbol_map[canonical_id]
        attr_dict.pop("gene_name", None)  # remove old GO-style 'gene_name'

        cols[8] = ";".join(f"{k}={v}" for k, v in attr_dict.items())
        gff_out.write("\t".join(cols) + "\n")

print(f"Annotated GFF written to: {OUTPUT_GFF}")