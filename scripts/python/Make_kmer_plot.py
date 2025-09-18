import sys
import matplotlib.pyplot as plt



dir="/Users/lizicheng/Desktop/Data/dog_genome/oyster/final_assembly/primary_polished_dedup_chr_qual/superscfd_merqury_report"
data="superscfd_chr_merqury_output.spectra-asm.hist"


kmer_counts = {}

with open(dir+"/"+data, 'r') as f:
    for line in f:
    # Skip header lines starting with '#'
        if line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) == 3:
            try:
                assembly_name = parts[0]
                count=int(parts[2])
                frequency=int(parts[1])
                if assembly_name not in kmer_counts:
                    kmer_counts[assembly_name] = {}
                kmer_counts[assembly_name][frequency] = count
            except ValueError:
                print(f"Skipping malformed line: {line.strip()}", file=sys.stderr)
        else:
            print(f"Skipping malformed line: {line.strip()}", file=sys.stderr)

frequencies = sorted(kmer_counts.keys())
counts = [kmer_counts[freq] for freq in frequencies]

    # Create the plot
plt.figure(figsize=(10, 6))  # Adjust figure size as needed

for assembly_name, data in kmer_counts.items():
    frequencies = sorted(data.keys())
    counts = [data[freq] for freq in frequencies]
    if assembly_name=="read-only":
    	label="WGS reads"
    else:
    	label="Assembly"
    plt.plot(frequencies, counts, marker='o', linestyle='-', markersize=4, label=label) # Added label


plt.title("21-mer Spectrum Plot")
plt.xlabel("21-mer multiplicity")
plt.ylabel("Count")
plt.xlim(xmin=0,xmax=80)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)  # Add gridlines
plt.legend()
plt.tight_layout() #Adjusts the padding
plt.show()

