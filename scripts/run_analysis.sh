#!/bin/bash
#SBATCH --comment=4164020200
#SBATCH --time=00-20:00:00
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=predict_EPI

#SBATCH --mail-type=END
#SBATCH --mail-user=pascal.duenk@wur.nl

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=12000
#SBATCH --exclude=fat001,fat002,fat100,fat101
#-----------------------------Environment, Operations and Job steps----

cd /lustre/nobackup/WUR/ABGC/duenk002/EPI_prediction

# on part of the data (chromosome 1-6)

# Rscript scripts/analysis.R "part" "log" "rand" "label"
# Rscript scripts/analysis.R "part" "original" "rand" "label"
# Rscript scripts/analysis.R "part" "all" "rand" "label"
# Rscript scripts/analysis.R "part" "log" "chromosome" "label"
# Rscript scripts/analysis.R "part" "original" "chromosome" "label"
# Rscript scripts/analysis.R "part" "all" "chromosome" "label"

# on the entire data 
# Rscript scripts/analysis.R "all" "log" "rand" "label"
# Rscript scripts/analysis.R "all" "original" "rand" "label"
# Rscript scripts/analysis.R "all" "all" "rand" "label"
# Rscript scripts/analysis.R "all" "log" "chromosome" "label"
# Rscript scripts/analysis.R "all" "original" "chromosome" "label"
# Rscript scripts/analysis.R "all" "all" "chromosome" "label"

Rscript scripts/analysis.R "all" "log" "promoter" "label"