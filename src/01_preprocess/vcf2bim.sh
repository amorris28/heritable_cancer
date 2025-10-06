#!/bin/bash
#SBATCH --account=p2521_tsd
#SBATCH --job-name=vcf2bim
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

# Process imputed vcf files.
# Make pvar files with all SNPs in the imputed vcf files.
# Make pgen/psam/pvar files applying filters on DR2 (imputation quality metrics) and MAF.


CHR=${SLURM_ARRAY_TASK_ID}
out_dir="data/interim/husk_bfiles"
vcf="/ess/p2521/data/durable/HuskMent/genotypes/imputed/v240619/vcf/husk_all.imputed.chr${CHR}.vcf.gz"
bfile="${out_dir}/husk_all.chr${CHR}"

mkdir -p ${out_dir}

#module load PLINK/2.00a2.3_x86_64
#plink2 --vcf ${vcf} dosage=DS \
#    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1024 error \
#    --make-bed --out ${bfile} --memory 60800 --threads 8

module load PLINK/1.90b_6.2-x86_64
plink --vcf ${vcf} --make-just-bim --out ${bfile}
