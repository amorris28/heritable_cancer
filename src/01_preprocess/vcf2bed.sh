#!/bin/bash
#SBATCH --account=p2521_tsd
#SBATCH --job-name=vcf2bed
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=7600M
#SBATCH --cpus-per-task=8
#SBATCH --array=3-22
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

# Process imputed vcf files.
# Make pvar files with all SNPs in the imputed vcf files.
# Make pgen/psam/pvar files applying filters on DR2 (imputation quality metrics) and MAF.

module load PLINK/2.00a2.3_x86_64

CHR=${SLURM_ARRAY_TASK_ID}
out_dir="data/interim/husk_bfiles"
vcf="/ess/p2521/data/durable/HuskMent/genotypes/imputed/v240619/vcf/husk_all.imputed.chr${CHR}.vcf.gz"
bfile="${out_dir}/husk_all.chr${CHR}"

mkdir -p ${out_dir}

plink2 --vcf ${vcf} dosage=DS \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1024 error \
    --make-bed --out ${bfile} --memory 60800 --threads 8
