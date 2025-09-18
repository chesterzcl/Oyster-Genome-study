#!/bin/bash
#SBATCH --job-name=depth_summary
#SBATCH --time=01:10:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# Load samtools module if needed
module load SAMtools  # Or adjust based on your cluster module

# Change to the directory with BAM files and position.txt
cd /home/zl436/palmer_scratch/bam/2nd_draft

# Output file
output_file="depth_at_selected_positions.tsv"

# Write header to output
echo -e "Sample\tPos_50102244\tPos_50102256\tPos_50102265" > "$output_file"

# Loop through all BAM files and extract depth at the 3 positions
for bam in *_filtered.bam; do
    echo -n "$bam" >> "$output_file"
    samtools depth -a -b position.txt "$bam" | awk '{printf "\t%s", $3} END {print ""}' >> "$output_file"
done
