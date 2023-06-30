mkdir 00_raw_data
mkdir 01_extract
mkdir 02_hisat2
mkdir 03_samtools
mkdir 04_counts

#### 6. umi-tools whitelist (cell barcode 8+6+6 in Read 1)
mv *gz ./00_raw_data
cd  00_raw_data
#rename "s/.merge.fq.gz/.fq.gz/" *
for i in ./*_1.fq.gz; do base=$(basename $i "_1.fq.gz");  /media/xionglab/data1/02_Software/umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTT){e<=1}.*'  --extract-method=regex --ed-above-threshold=correct --error-correct-threshold=2 --set-cell-number 1500 --knee-method=distance --log2stderr > ${base}_whitelist.txt --plot-prefix ${base}_expect_whitelist; done
for i in ./*_whitelist.txt; do base=$(basename $i "_whitelist.txt");  /media/xionglab/data1/02_Software/umi_tools extract --ignore-read-pair-suffixes --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCT){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>TTTTTTTTTT){e<=1}.*' --stdin ${base}_1.fq.gz  --stdout ../01_extract/${base}_1.extract.fq.gz  --read2-in  ${base}_2.fq.gz  --read2-out ../01_extract/${base}_2.extract.fq.gz   --error-correct-cell --extract-method=regex  --whitelist ${base}_whitelist.txt; done

cd ../01_extract
for i in `ls *_2.extract.fq.gz`
do base=$(basename $i "_2.extract.fq.gz")
cutadapt -a AAAAAAAA  -A AAAAAAAA -q 20 -O 8  --trim-n  -m 6  --max-n 0.1 -o ${base}_1.trimmed.fq.gz -p ${base}_2.trim.fq.gz ${base}_1.extract.fq.gz ${base}_2.extract.fq.gz
cutadapt -a CTGTCTCTTATACAC  -q 20 -O 8  --trim-n  -m 20  --max-n 0.1 -o ${base}_2.trimmed.fq.gz ${base}_2.trim.fq.gz
hisat2 -x /media/xionglab/data1/01_Database/03_Hisat2/mm10/genome -p 4 -U ${base}_2.trimmed.fq.gz -S ../02_hisat2/${base}.sam 2>../02_hisat2/${base}.mm10.align.log
samtools view -hbS -q 20 ../02_hisat2/${base}.sam | samtools sort -T ${base} - > ../03_samtools/${base}.bam
#featureCounts -a /media/xionglab/data1/01_Database/03_Hisat2/hisat2_gtf/refGene.mm10.gtf -o exon_assigned -R BAM ../03_samtools/${base}.bam -T 4 
featureCounts -a /media/xionglab/data1/01_Database/03_Hisat2/hisat2_gtf/mm10_Refseq.transcript.uniq2.gtf -t transcript -o ../03_samtools/${base}.gene_assigned -R BAM ../03_samtools/${base}.bam -T 6 
samtools sort ../03_samtools/${base}.bam.featureCounts.bam -o ../03_samtools/${base}.assigned_sorted.bam
samtools index ../03_samtools/${base}.assigned_sorted.bam
/media/xionglab/data1/02_Software/umi_tools count --wide-format-cell-counts --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I ../03_samtools/${base}.assigned_sorted.bam -S ../04_counts/${base}.counts.csv
done 