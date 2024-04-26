# Read count extraction from Zymo-Seq SwitchFree 3' mRNA library prep (Cat R3008) and NovaSeq X Plus 2 x 150bp seqeuncing

```
nano /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/fastqc_raw.sh
```

```  
#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae

module load all/FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cat /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/*.md5 > /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/Genohub_md5sum_list.txt

md5sum /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/*.gz > /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/URI_md5sum_list.txt

for file in /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/*.gz
do
fastqc $file --outdir /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/
done

multiqc /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/
```

```
sbatch /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/fastqc_raw.sh
```

scp hputnam@ssh3.hac.uri.edu://data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/multiqc_report.html /Users/hputnam/MyProjects/Poc_RAPID/RAnalysis/data/Tagseq_zymo_genohub/


## Setting up HPC folder for analysis

```
#set the main folder name
mkdir /data/putnamlab/hputnam/Pverr_Larvae_Devo
cd /data/putnamlab/hputnam/Pverr_Larvae_Devo
#make a raw data folder
mkdir raw
#move the data into the folder
cp /data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/*.gz raw/
#make folders for analysis
mkdir scripts
mkdir clean
mkdir refs

#Obtain the reference genome
#Reference Genome - https://academic.oup.com/gbe/article/12/10/1911/5898631
mkdir Pverr_Genome
cd Pverr_Genome/

#Genome scaffolds	fb4d03ba2a9016fabb284d10e513f873
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz

#Gene models (CDS)	019359071e3ab319cd10e8f08715fc71
wget http://pver.reefgenomics.org/download/Pver_genes_names_v1.0.fna.gz

#Gene models (proteins)	438f1d59b060144961d6a499de016f55
wget http://pver.reefgenomics.org/download/Pver_proteins_names_v1.0.faa.gz

#Gene models (GFF3)	614efffa87f6e8098b78490a5804c857
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.gff3.gz

#Full transcripts	76b5d8d405798d5ca7f6f8cc4b740eb2
wget http://pver.reefgenomics.org/download/Pver_transcriptome_v1.0.fasta.gz

md5sum *.gz > URI.download.md5


019359071e3ab319cd10e8f08715fc71  Pver_genes_names_v1.0.fna.gz
fb4d03ba2a9016fabb284d10e513f873  Pver_genome_assembly_v1.0.fasta.gz
614efffa87f6e8098b78490a5804c857  Pver_genome_assembly_v1.0.gff3.gz
438f1d59b060144961d6a499de016f55  Pver_proteins_names_v1.0.faa.gz
76b5d8d405798d5ca7f6f8cc4b740eb2  Pver_transcriptome_v1.0.fasta.gz

```

# Trimming reads

The intial QC showed that Read1 had high unique reads and Read2 had high amount of duplicated reads. This is from the library prep process where the polyA priming generates reads that are characterized as having redundancy by FastQC, whereas the UMIs on Read1 result in FastQC characterizing them as unique. 

The Zymo library prep kit for the Zymo-Seq SwitchFree 3' mRNA library prep (Cat R3008) recommends to only anlyze Read 2. I was worried about the amount of duplicated reads in Read2, so first I only analyzed Read1 which contains the P5 adapter, UMI + 6V, Oligo dT and then the insert. This requires both adapter trimming and trimming X number of bases in from the front of Read 1 to get to the Insert.


# Read 1 Analysis

```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R1.sh
```

```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/raw

#load modules
module load fastp/0.19.7-foss-2018b
module load all/FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim
array1=($(ls *R1_001.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean_${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
# fastqc the cleaned reads
        fastqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean_${i} --outdir /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/
done 

multiqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/

```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R1.sh
```

scp -r hputnam@ssh3.hac.uri.edu://data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/*.html /Users/hputnam/MyProjects/Poc_RAPID/RAnalysis/data/Tagseq_zymo_genohub/


## QC shows large T bias from bases 10 - 30 in read 1. I will run another trimming to remove the first 30bp of read 1

Also there is still adapter seqeunce present. Kit instructions say 

"Generally, it is sufficient to trim Read 2 with the Illumina® TruSeq® adapter
sequence: “AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT”."

In the unlikely event that Read 1 is required for downstream analysis,
please trim the random hexamer and the oligo dT sequences from the
beginning of Read 1 in addition to trimming the Illumina® TruSeq® adapter
sequence, “AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC”, from
the end of Read 1.

Here we will use AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R1_front30.sh
```

```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/

#load modules
module load fastp/0.19.7-foss-2018b
module load all/FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim
array1=($(ls *R1_001.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/mapped2/${i} --adapter_sequence=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim_front1 30 
# fastqc the cleaned reads
        fastqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/mapped2/${i} --outdir /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/
done 

multiqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30*

```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R1_front30.sh
```


# Test alignment Read 1

```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/mapping.sh
```

```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/

#load modules
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.9-foss-2018b


#hisat2-build -f /data/putnamlab/hputnam/Pverr_Larvae_Devo/refs/Pverr_Genome/Pver_genome_assembly_v1.0.fasta ./Pverr_Genome 
#echo "Referece genome indexed. Starting alingment" $(date)

	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_106_S51_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/106.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/106.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/106.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_107_S52_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/107.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/107.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/107.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_108_S53_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/108.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/108.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/108.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_109_S54_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/109.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/109.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/109.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_110_S55_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/110.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/110.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/110.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_111_S56_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/111.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/111.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/111.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_112_S57_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/112.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/112.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/112.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_113_S58_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/113.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/113.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/113.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_6_S47_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/6.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/6.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/6.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_7_S48_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/7.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/7.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/7.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_8_S49_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/8.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/8.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/8.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean/clean30_clean_9_S50_R1_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/9.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/9.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/9.sam


  
rm /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/*.sam


```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/mapping.sh
```

## Alignment results R1

|Sample.ID | overall alignment rate| aligned 0 times | aligned exactly 1 time | aligned >1 times | reads|
|---|---|---|---|---|---|
|106|79.79%|1114389 (20.21%)|3951617 (71.68%)|447170 (8.11%)|5513176 |
|107|72.84%|1553799 (27.16%)|3688108 (64.46%)|479373 (8.38%)|5721280 |
|108|66.91%|1331959 (33.09%)|2324323 (57.74%)|369070 (9.17%)|4025352 |
|109|79.90%|1368217 (20.10%)|4683859 (68.81%)|754715 (11.09%)|6806791|
|110|72.54%|1403566 (27.46%)|3119606 (61.04%)|587628 (11.50%)|5110800|
|111|79.42%|1197613 (20.58%|4043448 (69.50%)|576924 (9.92%)|5817985|
|112|79.67%|1199480 (20.33%)|4144181 (70.22%)|557717 (9.45%)|5901378|
|113|75.33%|1371601 (24.67%)|3746268 (67.38%)|442097 (7.95%)|5559966|
|6|80.93%|1111841 (19.07%)|4151500 (71.20%)|567392 (9.73%)|5830733|
|7|81.77%|1156119 (18.23%)|4467037 (70.43%)|719533 (11.34%)|6342689|
|8|84.11%|985584 (15.89%)|4498417 (72.55%)|716802 (11.56%)|6200803|
|9|84.80%|779842 (15.20%)|3799632 (74.05%)| 552017 (10.76%)|5131491|



# WORK IN PROGRESS


# Assembly with Stringtie 2


```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/assemble_transcripts.sh
```


```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped

#load modules
module load StringTie/2.1.4-GCC-9.3.0

array=($(ls *.bam)) #Make an array of bam files to assemble
 
for i in ${array[@]}; do 
	stringtie -p 8 --rf -e -G /data/putnamlab/hputnam/Pverr_Larvae_Devo/refs/Pverr_Genome/Pver_genome_assembly_v1.0.gff3 -o ${i}.gtf ${i}
done


```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/assemble_transcripts.sh
```


# Generate gene counts matrix
[prepDE.py](https://raw.githubusercontent.com/gpertea/stringtie/master/prepDE.py)


```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/prepDE_matrix.sh
```


```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped

#load modules
module load Python/2.7.15-foss-2018b 
module load StringTie/2.1.4-GCC-9.3.0
module load GffCompare/0.12.1-GCCcore-8.3.0


#make gtf list text file
ls *.bam.gtf > gtf_list.txt

stringtie --merge -p 8 -G /data/putnamlab/hputnam/Pverr_Larvae_Devo/refs/Pverr_Genome/Pver_genome_assembly_v1.0.gff3 -o ST_merged.gtf gtf_list.txt

gffcompare -r /data/putnamlab/hputnam/Pverr_Larvae_Devo/refs/Pverr_Genome/Pver_genome_assembly_v1.0.gff3 -G -o merged ST_merged.gtf 

array=($(ls *.bam))

for i in ${array[@]}; do 
	stringtie -p 8 --rf -e -G ST_merged.gtf -o ${i}.gtf ${i}
done


for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/prepDE.py -g Pverr_larvae_gene_count_matrix.csv -i listGTF.txt

```


```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/prepDE_matrix.sh
```

# The gene counts workflow is failing at the prepDE.py level with the following error

```
Problem parsing file /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped/106.bam.gtf at line:
:['Pver_Sc0000112_size552208', '.', 'transcript', '1269', '1387', '.', '+', '.', 'transcript_id "Pver_g6505.t1', '5_prime_partial true"; cov "0.0"; FPKM "0.000000"; TPM "0.000000";\n']
```


# Read 2 Data Analysis


```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R2.sh
```

```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/raw

#load modules
module load fastp/0.19.7-foss-2018b
module load all/FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim
array1=($(ls *R2_001.fastq.gz)) 

# fastp loop; trim the Read 2 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_${i} --adapter_sequence=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim_poly_x 6 -q 30 -y -Y 50 
# fastqc the cleaned reads
        fastqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_${i} --outdir /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/
done 

multiqc /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/

```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/trimQC_R2.sh
```


scp hputnam@ssh3.hac.uri.edu://data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/multiqc_report.html /Users/hputnam/MyProjects/Poc_RAPID/RAnalysis/data/Tagseq_zymo_genohub/




# Test Mapping Read 2

```
nano /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/R2_mapping.sh
```

```
#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=400GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/hputnam/Pverr_Larvae_Devo/

#load modules
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.9-foss-2018b

#hisat2-build -f /data/putnamlab/hputnam/Pverr_Larvae_Devo/refs/Pverr_Genome/Pver_genome_assembly_v1.0.fasta ./Pverr_Genome 

	#hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_106_S51_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/106.sam
        #samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/106.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/106.sam

hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_107_S52_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/107.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/107.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/107.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_108_S53_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/108.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/108.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/108.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_109_S54_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/109.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/109.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/109.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_110_S55_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/110.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/110.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/110.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_111_S56_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/111.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/111.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/111.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_112_S57_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/112.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/112.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/112.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_113_S58_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/113.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/113.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/113.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_6_S47_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/6.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/6.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/6.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_7_S48_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/7.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/7.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/7.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_8_S49_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/8.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/8.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/8.sam


	hisat2 --rna-strandness F -p 8 --dta -x Pverr_Genome -U /data/putnamlab/hputnam/Pverr_Larvae_Devo/clean2/clean_9_S50_R2_001.fastq.gz -S /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/9.sam
        samtools sort -@ 8 -o /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/9.bam /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/9.sam


  
rm /data/putnamlab/hputnam/Pverr_Larvae_Devo/mapped2/*.sam

```

```
sbatch /data/putnamlab/hputnam/Pverr_Larvae_Devo/scripts/R2_mapping.sh
```

## Alignment results Read2

|Sample.ID | overall alignment rate| aligned 0 times | aligned exactly 1 time | aligned >1 times | reads|
|---|---|---|---|---|---|
|106|68.97%|2617559 (31.03%)|5255055 (62.30%)|562400 (6.67%)|8435014 |
|107|58.23%|3818475 (41.77%)|4566278 (49.95%)|756429 (8.27%)|9141182|
|108|39.12%|4575629 (60.88%)|2443027 (32.51%)|496571 (6.61%)|7515227|
|109|67.74%|3190933 (32.26%)|5902975 (59.68%)|797413 (8.06%)|9891321|
|110|51.71%|4423586 (48.29%)|3758092 (41.03%)|978454 (10.68%)|9160132|
|111|69.42%|2660851 (30.58%)|5362514 (61.64%)|676534 (7.78%)|8699899|
|112|69.42%|2679906 (30.58%)|5476030 (62.49%)|607125 (6.93%)|8763061|
|113|64.29%|3033008 (35.71%)|4893828 (57.62%)|565831 (6.66%)|8492667|
|6|54.29%|4173872 (45.71%)|4521417 (49.51%)|436760 (4.78%)|9132049|
|7|54.73%|4307010 (45.27%)|4615189 (48.51%)|591322 (6.22%)|9513521|
|8|72.38%|2567051 (27.62%)|5914595 (63.65%)|811121 (8.73%)|9292767|
|9|75.36%|1902138 (24.64%)|5111941 (66.21%)| 707090 (9.16%)|7721169|

