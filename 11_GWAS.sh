#!/bin/bash
#! What is this job called?
#SBATCH --job-name=black_GWAS
#! What should I call the reports
#SBATCH --output="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/GWAS_%A"
#SBATCH --error="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/GWAS_%A"
#! Which project should be charged:
#SBATCH -A BURGESS-SL3-CPU
#! Which partition/cluster am I using?
#SBATCH -p skylake
#! How many nodes should be allocated? 
#SBATCH --nodes=1
#! How many tasks will there be in total? By default SLURM
#SBATCH --ntasks=1
#! How much memory in MB is required _per node_? 
##SBATCH --mem=5Gb
#! How many cpus per task
#SBATCH --cpus-per-task=23
#! How much wallclock time will be required?
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL

module load gcc/5
module load R/4.2.0
module load pandoc/2.0.6

Rscript "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/GWAS_blackcohort.R" --exposure="ischtia" --exposure_family="logistic" 
