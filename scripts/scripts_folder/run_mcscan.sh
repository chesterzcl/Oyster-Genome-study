#!/bin/bash

# Load your module if needed
module load GCC

# Optional: move to your working directory
cd /home/zl436/palmer_scratch/script/mcscan_dir

input_name=MG_CV_final/MG_vs_CV_final
# Run MCScanX
MCScanX/MCScanX ${input_name}
