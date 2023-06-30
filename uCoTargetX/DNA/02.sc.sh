######## single-cell 
#手动把02_splitT5_data里的文件按mark合并，拖到01_splitT7_data中
cd /media/helab/data1/haiqing/20230202_QW_293AD/merge/plan1/DNA
mkdir 04_extract
mkdir 05_bowtie2
mkdir 06_sort_uni_bam
mkdir 07_rmdup_bam
mkdir -p 08_split/01_H3K27ac/sc
mkdir -p 08_split/03_H3K27me3/sc

#### 6. umi-tools whitelist (cell barcode 8+6+6 in Read 1)
cd ./01_splitT7_data
for i in ./*.1.fq.gz; do base=$(basename $i ".1.fq.gz");  umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG{e<=1}).*'  --extract-method=regex --ed-above-threshold=correct  --error-correct-threshold=2 --set-cell-number 1000 --log2stderr --knee-method=distance > ${base}_whitelist.txt --plot-prefix ${base} ; done
cd ..

#### 7. umi-tools extract
cd 01_splitT7_data
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  umi_tools extract --ignore-read-pair-suffixes --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG{e<=1}).*' --stdin ${base}.1.fq.gz  --stdout ../04_extract/${base}_1.extract.fq.gz  --read2-in  ${base}.2.fq.gz  --read2-out ../04_extract/${base}_2.extract.fq.gz  --error-correct-cell --extract-method=regex  --whitelist ${base}_whitelist.txt; done
cd ../04_extract/

#### 8. cutadapt and mapping & 
# barcode info was recorded in read 1 name 
# read 1 has been trimmed by umi-tools
for i in `ls *_2.extract.fq.gz`
do base=$(basename $i "_2.extract.fq.gz")
cutadapt -u 52 -o ${base}.cut.2.extract.fq.gz ${base}_2.extract.fq.gz
cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -B CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 10 --max-n 0.1 --trim-n -o ${base}.trimmed.1.fq.gz -p ${base}.trimmed.2.fq.gz ${base}_1.extract.fq.gz ${base}.cut.2.extract.fq.gz
bowtie2 -x /media/helab/data1/00_public/database/genomes/hg19_bowtie2/hg19 -1 ./${base}.trimmed.1.fq.gz -2 ./${base}.trimmed.2.fq.gz --very-sensitive-local --no-unal -X 3000 -S ../05_bowtie2/${base}.hg19.sam -p 8 2>../05_bowtie2/${base}.align.log
samtools view -hbS -q 10 ../05_bowtie2/${base}.hg19.sam | samtools sort -T -> ../06_sort_uni_bam/${base}.hg19.unique.bam
java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=../06_sort_uni_bam/${base}.hg19.unique.bam o=../07_rmdup_bam/${base}.hg19.rmdup.bam M=../07_rmdup_bam/${base}.hg19.txt
samtools view ../07_rmdup_bam/${base}.hg19.rmdup.bam | sed 's/_/\tCB:Z:/'  - | sed 's/_/\tTB:Z:/' - |awk '{print$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$2"\t"$3}' -> ../07_rmdup_bam/${base}_CB_TB.bam
samtools view ../07_rmdup_bam/${base}.hg19.rmdup.bam -H | cat - ../07_rmdup_bam/${base}_CB_TB.bam | samtools sort -t CB - -o ../07_rmdup_bam/${base}_CB_TB_sorted_tags.bam
done 

#### 9. single-cell processing
cp `find ../07_rmdup_bam -name "*_CB_TB_sorted_tags.bam" |grep "H3K27ac"` ../08_split/01_H3K27ac
cp `find ../07_rmdup_bam -name "*_CB_TB_sorted_tags.bam" |grep "H3K27me3"` ../08_split/03_H3K27me3

python ../split_script_v2_k27ac.py &
python ../split_script_v2_k27.py &
