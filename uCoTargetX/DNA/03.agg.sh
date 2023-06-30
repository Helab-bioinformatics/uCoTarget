mkdir -p 03_aggregate/00_raw_data
mkdir -p 03_aggregate/01_cutadapt
mkdir -p 03_aggregate/02_bowtie2
mkdir -p 03_aggregate/03_sort_uni_bam
mkdir -p 03_aggregate/04_rmdup_bam
mkdir -p 03_aggregate/05_bw

#注意，run这一步需要保证01_splitT7_data子文件里没有fq.gz,否则会重复merge
cat `find ./01_splitT7_data -name "*.1.fq.gz" |grep "H3K27ac"|sort` > 03_aggregate/00_raw_data/merge.H3K27ac.1.fq.gz
cat `find ./01_splitT7_data -name "*.2.fq.gz" |grep "H3K27ac"|sort` > 03_aggregate/00_raw_data/merge.H3K27ac.2.fq.gz

cat `find ./01_splitT7_data -name "*.1.fq.gz" |grep "H3K27me3"|sort` > 03_aggregate/00_raw_data/merge.H3K27me3.1.fq.gz
cat `find ./01_splitT7_data -name "*.2.fq.gz" |grep "H3K27me3"|sort` > 03_aggregate/00_raw_data/merge.H3K27me3.2.fq.gz

cd 03_aggregate/00_raw_data
for i in *gz 
do
fastqc $i 
done&

for i in `ls *.2.fq.gz`
do base=$(basename $i ".2.fq.gz")
cutadapt -u 98 -o ../01_cutadapt/${base}.cut.1.fq.gz ${base}.1.fq.gz
cutadapt -u 52 -o ../01_cutadapt/${base}.cut.2.fq.gz ${base}.2.fq.gz
cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 30 --max-n 0.1 --trim-n -o ../01_cutadapt/${base}.trimmed.1.fq.gz -p ../01_cutadapt/${base}.trimmed.2.fq.gz ../01_cutadapt/${base}.cut.1.fq.gz ../01_cutadapt/${base}.cut.2.fq.gz
bowtie2 -x /media/helab/data1/00_public/database/genomes/hg19_bowtie2/hg19 -1 ../01_cutadapt/${base}.trimmed.1.fq.gz -2 ../01_cutadapt/${base}.trimmed.2.fq.gz --very-sensitive-local --no-unal -X 3000 -S ../02_bowtie2/${base}.hg19.sam -p 8 2>../02_bowtie2/${base}.align.log
samtools view -hbS -q 20 ../02_bowtie2/${base}.hg19.sam | samtools sort -T -> ../03_sort_uni_bam/${base}.hg19.unique.bam
java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=../03_sort_uni_bam/${base}.hg19.unique.bam o=../04_rmdup_bam/${base}.hg19.rmdup.bam M=../04_rmdup_bam/${base}.hg19.txt
samtools index ../04_rmdup_bam/${base}.hg19.rmdup.bam
done  


cd ../04_rmdup_bam
for i in ./*.rmdup.bam
do
base=$(basename $i ".rmdup.bam")
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i  -e 300  --smoothLength 500 -o  ../05_bw/${base}.ext300.smo500.bw -p 8
done

cd ../03_sort_uni_bam
for x in ./*.unique.bam
do 
base=$(basename $x ".unique.bam")
#samtools index $x
samtools view -c  ${base}.unique.bam  >./${base}.uni.totolreads.log
done


cd ../04_rmdup_bam
for j in ./*.rmdup.bam
do 
base=$(basename $j ".rmdup.bam")
samtools view -c ${base}.rmdup.bam  >./${base}.totolreads.log
done
