
#此版本仅用于测试，正式版本会在不久后更新

#准备数据 放在seq文件夹下面

#groEL     TCCGATTACGAYCGYGAGAAGCT   CSGCYTCGGTSGTCAGGAACAG


require cutadapt
require vsearch



#对于单端数据
time for i in `tail -n+2 metadata.txt|cut -f1`;do
 cutadapt -e 0.1 -g TCCGATTACGAYCGYGAGAAGCT seq/${i}.fq >  trum1/${i}.fq &
done 

time for i in `tail -n+2 metadata_lu.txt|cut -f1`;do
 cutadapt -e 0.1 -a CTGTTCCTGACSACCGADGCSG trum1/${i}.fq >  trum2/${i}.fq &
done 

time for i in `tail -n+2 metadata_lu.txt|cut -f1`;do
 usearch11 -fastx_relabel trum2/${i}.fq -fastqout temp/${i}.fq -prefix ${i}.
done &




#对于双端数据
#切除引物
mkdir trum
time for i in `tail -n+1 metadata.txt|cut -f1`;do
cutadapt -u -10 -e 0.1 -g TCCGATTACGAYCGYGAGAAGCT seq/${i}.R1.fq.gz >  trum/${i}_R1.fq.gz &
cutadapt -u -10 -e 0.1 -g CSGCYTCGGTSGTCAGGAACAG seq/${i}.R2.fq.gz >  trum/${i}_R2.fq.gz &
 done 

# # 批量改名，需要有单端fastq文件，且解压(usearch不支持压缩格式)
mkdir temp
time for i in `tail -n+1 metadata.txt|cut -f1`;do
 vsearch --fastq_mergepairs trum/${i}_R1.fq.gz --reverse trum/${i}_R2.fq.gz \
 --fastqout temp/${i}.fq --relabel ${i}. &
done 


 #合并所有样品至同一文件
 cat temp/*.fq > temp/all.fq



#查看序列质量情况
seqkit stats temp/all.fq


# 质控
nohup vsearch --fastx_filter temp/all.fq \
 --fastq_maxee_rate 0.01 \
 --fastq_minlen 350 \
 --fastq_qmax 93 \
 --fastaout temp/filtered.fa &




#去冗余
nohup vsearch --derep_fulllength temp/filtered.fa \
 --minuniquesize 1 --sizeout --relabel Uni_ \
 --output temp/uniques.fa  &



#去嵌合体，denovo开销较大，如果文件较大，则不要考虑该步骤。
vsearch --uchime_deno temp/uniques.fa \
--sizein --sizeout --fasta_width 0 \
    --threads 20 \
   --nonchimeras temp/denovo.nonchimeras.fa

#去嵌合体，有参模板
vsearch --uchime_ref temp/denovo.nonchimeras.fa    \
  --db ./db/groEL_db.fa       \
  --thread 72 \
   --sizein --sizeout --fasta_width 0   \
   --nonchimeras temp/ref.nonchimeras.fa 



#聚类ASV
vsearch --cluster_unoise temp/ref.nonchimeras.fa \
--minsize 4 \
--centroids temp/zotus.fa \
--relabel ASV_ 




mkdir result
#修改序列名：Zotu为改为ASV方便识别
sed 's/Zotu/ASV_/g' temp/zotus.fa > result/otus.fa
head -n 2 result/otus.fa


#生成otu表格
vsearch --usearch_global temp/filtered.fa \
      --db result/otus.fa \
      --id 0.97 --threads 12 \
    	--otutabout result/otutab.txt 



#物种注释
./tools/taxon --db ./db/groEL_db.fa --infile result/otus.fa --outfile result.txt


