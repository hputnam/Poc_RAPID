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

scp hputnam@ssh3.hac.uri.edu://data/putnamlab/KITT/hputnam/20240424_Pverrucosa_larvae/multiqc_report.html Desktop/