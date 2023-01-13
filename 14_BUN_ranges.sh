#!/bin/bash
#SBATCH --time=12:0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --job-name=qctool2_black
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
#chr=${SLURM_ARRAY_TASK_ID}


# make new bgen of snps

#from scratch


qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr1.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 01:92559-94559 -incl-range 01:62181161-62183161 -incl-range 01:116841952-116843952 -incl-range 01:155380399-155382399 -incl-range 01:166431825-166433825  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges1.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr2_#.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 02:13455372-13457372 -incl-range 02:126648963-126650963 -incl-range 02:139139815-139141815 -incl-range 02:205900027-205902027 -incl-range 02:214697317-214699317 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges2.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr3.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 03:106040996-106042996 -incl-range 03:113187304-113189304  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges3.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr4.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 04:29031872-29033872 -incl-range 04:91034013-91036013  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges4.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr5.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 05:53478171-53480171 -incl-range 05:54401043-54403043 -incl-range 05:94288773-94290773 -incl-range 05:108688489-108690489 -incl-range 05:118537576-118539576 -incl-range 05:157057055-157059055 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges5.bgen" -ofiletype bgen_v1.2 -bgen-bits 8



qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr7.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 07:39069119-39071119 -incl-range 07:111430123-111432123 -incl-range 07:146417975-146419975   -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges7.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr8.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 08:37795973-37797973 -incl-range 08:117690075-117692075 -incl-range 08:123031487-123033487 -incl-range 08:141193073-141195073  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges8.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr10.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 10:4780174-4782174 -incl-range   10:106576220-106578220 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges10.bgen" -ofiletype bgen_v1.2 -bgen-bits 8



qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr11.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 11:24940022-24942022 -incl-range  11:39606225-39608225 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges11.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr13.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 13:45764379-45766379 -incl-range 13:72489638-72491638 -incl-range 13:85281031-85283031 -incl-range  13:105340518-105342518  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges13.bgen" -ofiletype bgen_v1.2 -bgen-bits 8



qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr14.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 14:44910041-44912041 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges14.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr15.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 15:41626530-41628530 15:80808194-80810194   -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges15.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr16.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 16:1002244-1004244 -incl-range 16:20577280-20579280 -incl-range 16:20577281-20579281 -incl-range 16:79064987-79066987  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges16.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr17.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 17:5743718-5745718  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges17.bgen" -ofiletype bgen_v1.2 -bgen-bits 8


qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr18.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 18:59635815-59637815 -incl-range 18:76130881-76132881 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges18.bgen" -ofiletype bgen_v1.2 -bgen-bits 8


qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr19.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 19:53368270-53370270 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges19.bgen" -ofiletype bgen_v1.2 -bgen-bits 8


qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr20.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -incl-range 20:1745883-1747883 -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges20.bgen" -ofiletype bgen_v1.2 -bgen-bits 8

# combine

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/BUN_keyranges#.bgen" -s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample"  -og "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/BUN_keyranges_all.bgen" -ofiletype bgen_v1.2 -bgen-bits 8


# make stats 

qctool -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/BUN_keyranges_all.bgen"-s "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample" -snp-stats -osnp "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/data/BUN_keyranges_snpstats.txt" 


# make index

module load ceuadmin/bgen

bgenix -g "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/BUN_keyranges_all.bgen" -index -clobber




