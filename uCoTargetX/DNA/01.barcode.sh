#### 1. create 00_raw_data folder for raw .fq.gz data
# fill the sample_info excel before analysis

#### 2. mkdir 
cd ./00_raw_data
rename "s/.merge.fq.gz/.fq.gz/" *
mkdir ../01_splitT7_data
mkdir ../02_splitT5_data
mkdir ../03_aggregate

#### 3. split Ab(T7) barcode
Rscript  ../create.T7.T5.R #revise the PATH of working directory in R
sh ../create.splitT7.sh
cd ../01_splitT7_data
rm ./*_noid.2.fq.gz
sh ../create.pickT7.sh

#### 4. split T5 barcode
sh ../create.splitT5.sh
rm ../02_splitT5_data/*_noid.1.fq.gz
cd ../02_splitT5_data
sh ../create.pickT5.sh

######## 5. aggregate analysis
# revise group name 
cd ..
mkdir -p 03_aggregate/00_raw_data
mkdir -p 03_aggregate/01_cutadapt
mkdir -p 03_aggregate/02_bowtie2
mkdir -p 03_aggregate/03_sort_uni_bam
mkdir -p 03_aggregate/04_rmdup_bam
mkdir -p 03_aggregate/05_bw