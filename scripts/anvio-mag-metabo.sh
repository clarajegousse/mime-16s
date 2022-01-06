salloc -N 1
ssh $SLURM_NODELIST
echo $HOSTNAME

# to insure work with python3
source /users/home/cat3/.bashrc
conda activate anvio-master

WD=/users/home/cat3/projects/mime-metabo
mkdir -p $WD
cd $WD

# ---- start with Synechococcus -----

mkdir METABAT__194-Synechococcus
cd METABAT__194-Synechococcus

wget 'ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737978/METABAT__194-surface-contigs.fa.gz'
gunzip METABAT__194-surface-contigs.fa.gz

anvi-gen-contigs-database -f METABAT__194-surface-contigs.fa -o METABAT__194-contigs.db -n Synechococcus -T 10 

anvi-run-kegg-kofams -c METABAT__194-contigs.db -T 10

anvi-estimate-metabolism -c METABAT__194-contigs.db --only-complete
