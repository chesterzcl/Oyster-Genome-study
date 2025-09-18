#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --time=05:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

module load Python

source ~/my_python_module/bin/activate

python3 summarize_ld.py /home/zl436/palmer_scratch/script/ld_analysis/ld_500kb_r2.ld
