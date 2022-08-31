#!/bin/bash
#SBATCH --job-name=timing_benchmark
#SBATCH --output=timing_benchmark_%j
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8

... # activate SECOND env

# Execution
python -u timing_mp.py
