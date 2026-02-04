#!/bin/bash
#BSUB -J GWAS1
#BSUB -n 4
#BSUB -R rusage[mem=20G]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q high

cd /public/home/zhaoli/network-driven-gwas/kernel/gwas-new
mkfifo fd1;                          #################
exec 4<>fd1;                         #################
thred=4;
{
        for ((i=1;i<=${thred};i++));do
                echo;
        done;
}>&4                                 #################

/public/home/software/opt/bio/software/TASSEL/tassel5/run_pipeline.pl -Xmx20g \
  -fork1 -h /public/home/zhaoli/network-driven-gwas/kernel/gwas-new/merge_all_filter_sorted-453-1.hmp.txt \
  -filterAlign -filterAlignMinFreq 0 \
  -fork2 -r /public/home/zhaoli/network-driven-gwas/kernel/gwas-new/100grainweight.txt \
  -fork3 -q /public/home/zhaoli/network-driven-gwas/plant-high/gwas/509-pop-structure.txt \
  -fork4 -k /public/home/zhaoli/network-driven-gwas/plant-high/gwas/509-kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 \
  -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export output-1/"$a" -runfork1 -runfork2 -runfork3 -runfork4 
  
rm -f fd1; 
