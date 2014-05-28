for i in `cat filenames.txt`
do
echo "bash ~/ATAC_seq/mapping.sh  ${i} /mnt/lustre/home/nbanovich/bowtie_indicies/hg19" | qsub -l h_vmem=14g  -wd `pwd`
done
