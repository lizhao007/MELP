#BSUB -J emmax-up
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q high

#plink --vcf merge_all_filter_sorted-174.vcf --indep-pairwise 100 50 0.2 --out snp.int0.8maf0.05 --allow-extra-chr --make-bed #过滤SNP数据窗口大小为100kb, 每次扫描 50个SNP,过滤标准为LD 0.2
#plink --bfile snp.int0.8maf0.05 --extract snp.int0.8maf0.05.prune.in --out prunData --recode 12 --allow-extra-chr #admixture格式，生成群体结构文件
#plink --bfile snp.int0.8maf0.05 --recode 12 --output-missing-genotype 0 --transpose --out snp_miss0.1_DP20 #生成emmax格式文件
#./emmax-kin-intel64 snp_miss0.1_DP20 -v -d 10 #利用emmax生成亲缘关系矩阵（kinship）
#plink --vcf merge_all_filter_sorted-174.vcf --recode 12 transpose --out emmax_in --allow-extra-chr --chr-set 10

#174个材料all-gene-2k1k内SNP
#/public/home/software/opt/bio/software/TASSEL/tassel5/run_pipeline.pl -Xmx10g -fork1 -h all-gene-174.hmp.txt -export -exportType VCF -runfork1 
#module load VCFtools/0.1.16
###这行命令转换的格式有错，tped文件中包含字母，而不是012，换个方式
#module load VCFtools/0.1.16
##vcftools --vcf all-gene-174.vcf --plink --out snp.int0.8maf0.05-600w.plk#module load plink/1.9
#plink --file snp.int0.8maf0.05-600w.plk --recode transpose --out snp.int0.8maf0.05-600w
####问题解决
#module load plink/1.9
#plink --vcf all-gene-174.vcf --indep-pairwise 100 50 0.2 --allow-extra-chr --make-bed --out snp.int0.8maf0.05-600w.plk
#plink --bfile snp.int0.8maf0.05-600w.plk --recode 12 --output-missing-genotype 0 --transpose --out snp.int0.8maf0.05-600w.plk
#for gene in $(cat gene.id); do ./emmax-intel64 -v -d 10 -t /public/home/zhaoli/network-driven-gwas/plant-high/snp.int0.8maf0.05-600w.plk -p ./phenotype/phenotype_${gene}.txt -k /public/home/zhaoli/network-driven-gwas/plant-high/snp_miss0.1_DP20.aBN.kinf -c /public/home/zhaoli/network-driven-gwas/plant-high/174-pop-structure.txt -o ./emmax-output/plant-high-${gene}.txt;done

#./emmax-intel64 -v -d 10 -t /public/home/zhaoli/network-driven-gwas/plant-high/snp.int0.8maf0.05-600w.plk -p ./phenotype/phenotype_Zm00001d054079.txt -o ./emmax-output/plant-high-Zm00001d054079.txt
#down 
for gene in $(cat gene.id); do ./emmax-intel64 -v -d 10 -t /public/home/zhaoli/perennial/suanfa/gwas/2022-gb-zeamap/index/1000w-181.plk -p ./phenotype-rna/phenotype_${gene}.txt -k /public/home/zhaoli/perennial/suanfa/gwas/2022-gb-zeamap/kinship/1000w-181.aBN.kinf -c /public/home/zhaoli/perennial/suanfa/gwas/2022-gb-zeamap/structure/181-structure-k9.txt -o ./emmax-output-rna/kernel-${gene}.txt;done


