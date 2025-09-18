import pandas as pd
import os

os.chdir("/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly_v2")


# ---- 1) Load eggNOG annotation ----
# adjust sep and column names based on your file
eggnog = pd.read_csv("eggnog_annotation.emapper.annotations", sep="\t")

# Map: query ID -> description
desc_map = dict(zip(eggnog["query"], eggnog["description"]))

# ---- 2) Process GFF ----
with open("braker.gff", "r") as infile, open("braker_eggnog_annotated.gff", "w") as outfile:
    for line in infile:
        if line.startswith("#"):   # keep headers unchanged
            outfile.write(line)
            continue
        
        parts = line.strip().split("\t")
        if len(parts) < 9:
            outfile.write(line)
            continue
        
        attributes = parts[8]
        
        # Extract ID (some GFFs use ID=gene or ID=transcript)
        gene_id = None
        for field in attributes.split(";"):
            if field.startswith("ID="):
                gene_id = field.replace("ID=", "")
                break
        
        # If we have annotation for this gene_id, add it
        if gene_id and gene_id in desc_map:
            desc = desc_map[gene_id].replace(";", ",")  # avoid breaking GFF format
            attributes += f";eggNOG_description={desc}"
        
        parts[8] = attributes
        outfile.write("\t".join(parts) + "\n")