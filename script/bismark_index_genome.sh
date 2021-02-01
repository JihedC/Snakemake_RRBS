#! /bin/bash

#SBATCH  --job-name=plotprofile
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.chouaref@lumc.nl
#SBATCH -t 48:00:00
#SBATCH --mem=60000

echo START job at `date`
bismark_genome_preparation --verbose .

echo DONE at `date `
