#!/bin/bash
#SBATCH --time=12:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=index2
#FILENAME: qctool_cutout
#AUTHOR : Amy Mason
#PURPOSE: create subset of the UKbiobank bgen files to enable running snptest (single output file)
#INPUT: input_snps, a space seperated text file listing snps to extract. They must be in the form chr:pos with leading zeros if chr<10 ( eg 06:160493099)
#       output_name, what to call the output files
#       output_dir, where to put the output files (this will be set to the bgen folder in your project file)
#OUTPUT: one .bgen file, containing the extracted snps for that chr in .bgen format


module load gcc/5
module load qctool/v2.0.5

# choose chr number
chr=${SLURM_ARRAY_TASK_ID}

# make index

module load ceuadmin/bgenix/1.0.2

bgenix -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/black_MAFf_chr${chr}.bgen" -index