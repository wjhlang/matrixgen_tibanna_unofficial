#!/bin/bash
#SBATCH --job-name=all_matrix
#SBATCH --mail-user=wjhlang@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=remills1
#SBATCH --partition=standard
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH -o %x_%j.out

chroms=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22')

# Iterate the chroms using for loop
for chrom in ${chroms[@]}; do
 filename="${chrom}.sh"
 content1="#!/bin/bash
#SBATCH --job-name=run_${chrom}
#SBATCH --mail-user=wjhlang@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --account=remills1
#SBATCH --partition=standard
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH -o %x_%j.out"
 content2="source /home/wjhlang/miniconda3/etc/profile.d/conda.sh
conda activate umich
python -u scripts/main.py ${chrom}"

 echo "$content1" >> "$filename"
 echo "$content2" >> "$filename"
 sbatch "$filename"
done


