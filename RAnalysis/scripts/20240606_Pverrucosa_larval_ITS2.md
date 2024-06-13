# Pverrucosa larval ITS2


[Huffmyer Symportal NOtebook](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-ITS2-Analysis-Part-2/)

Symportal Key
vM2MzoSWwuxc2nOP0zOCM7ixG73h8eQeIggWlDV8ENTH8BMgnu(base)


```
interactive
module load Miniconda3/4.9.2
module load SymPortal

conda env create -f $EBROOTSYMPORTAL/symportal_env.yml 

```

# Create a reference database

```
#!/bin/bash
#SBATCH --job-name="SP_reference"
#SBATCH -t 500:00:00
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ashuffmyer@uri.edu
#SBATCH -D /data/putnamlab/ashuffmyer/SymPortal
#SBATCH --exclusive

module load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"
conda activate symportal_env

module load SymPortal/0.3.21-foss-2020b
module unload

export PYTHONPATH=/data/putnamlab/hputnam/Pverr_ITS2/SymPortal/:/data/putnamlab/hputnam/Pverr_ITS2/SymPortal/lib/python3.7/site-packages:$PYTHONPATH

export PATH=/data/putnamlab/hputnam/Pverr_ITS2/SymPortal/:/data/putnamlab/hputnam/Pverr_ITS2/SymPortal/bin:$PATH

python3 manage.py migrate

python3 populate_db_ref_seqs.py

module unload SymPortal/0.3.21-foss-2020b

echo "Mission Complete!" $(date)
```