#!/bin/bash
#SBATCH --job-name=braker
#SBATCH --time=1:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL

apptainer pull braker3.sif docker://teambraker/braker3:latest

mkdir -p $PWD/augustus_config

apptainer exec --bind $PWD:/data braker3.sif bash -c 'mkdir -p /data/augustus_config && cp -r /opt/Augustus/config/* /data/augustus_config/'

