#!/bin/bash
#SBATCH --mem=262144
#SBATCH --nodelist=moab
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nl91@duke.edu

date
hostname

cd /luscinia/nl91/kwave/scratch/pwr2p0
python kwave_fft.py
