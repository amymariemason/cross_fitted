#!/bin/bash
#! What is this job called?
#SBATCH --job-name=instrument
#! What should I call the reports
##SBATCH --output="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/instrument_%A_%a"
##SBATCH --error="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/instrument_%A_%a"
#! Which project should be charged:
#SBATCH -A BURGESS-SL3-CPU
#! Which partition/cluster am I using?
#SBATCH -p skylake
#! How many nodes should be allocated? 
#SBATCH --nodes=1
#! How many tasks will there be in total? By default SLURM
#SBATCH --ntasks=1
#! How much memory in MB is required _per node_? 
##SBATCH --mem=10Gb
#! How many cpus per task
#SBATCH --cpus-per-task=23
#! How much wallclock time will be required?
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL




module load gcc/5
module load R/4.2.0


Rscript "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/Mr_fitting.R" message_warning.R >& /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/step_10_output.txt


