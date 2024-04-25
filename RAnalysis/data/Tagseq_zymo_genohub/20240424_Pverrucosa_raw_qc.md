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

The intial QC showed that the read information is inverted from what is described in the Zymo library prep kit for the Zymo-Seq SwitchFree 3' mRNA library prep (Cat R3008).  So here we are analyzing Read 1 which contains the data and not Read 2 which contains the UMI and mostly adapter. 


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

scp -r hputnam@ssh3.hac.uri.edu://data/putnamlab/hputnam/Pverr_Larvae_Devo/*.html /Users/hputnam/MyProjects/Poc_RAPID/RAnalysis/data/Tagseq_zymo_genohub/


/data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/*_adapterTrimmed.fastq.gz