salloc -N 1
ssh $SLURM_NODELIST
echo $HOSTNAME

# to insure work with python3
source /users/home/cat3/.bashrc
conda activate anvio-master

WD=/users/home/cat3/projects/mime-mg
mkdir -p $WD
cd $WD

# ----- download TARA metagenomic datasets -----

# /users/home/cat3/projects/mime-mg/data/tara
mkdir -p $WD/data/tara
cd $WD/data/tara
#
# echo -e "ERR3589592\nERR3589586\nERR3589582" > acc.txt
#
# cat $WD/data/tara/acc.txt | while read -r acc
# do
# echo $acc
# l=$(expr length $acc)
#
# if [[ "$l" == 10 ]]; then
# wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/00'${acc:9:10}'/'$acc'/'$acc'_1.fastq.gz'
# wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/00'${acc:9:10}'/'$acc'/'$acc'_2.fastq.gz'
# else
# wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/'$acc'/'$acc'_1.fastq.gz'
# wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/'$acc'/'$acc'_2.fastq.gz'
# fi
# done
#
# echo -e "sample\tr1\tr2
# TARA_158\tERR3589592_1.fastq.gz\tERR3589592_2.fastq.gz
# TARA_163\tERR3589586_1.fastq.gz\tERR3589586_2.fastq.gz
# TARA_210\tERR3589582_1.fastq.gz\tERR3589582_2.fastq.gz" > samples.txt
#
# iu-gen-configs samples.txt -o $WD/data/tara
# iu-filter-quality-minoche $WD/data/tara/TARA_158.ini --ignore-deflines
# iu-filter-quality-minoche $WD/data/tara/TARA_163.ini --ignore-deflines
# iu-filter-quality-minoche $WD/data/tara/TARA_210.ini --ignore-deflines

sample="TARA_158"
R1=$sample'-QUALITY_PASSED_R1.fastq.gz'
R2=$sample'-QUALITY_PASSED_R1.fastq.gz'
flash $R1 $R2 -x 0.05 -m 20 -M 150 -z -o '$sample' 2>&1 | tee $sample'-flash.log'

sample="TARA_163"
R1=$sample'-QUALITY_PASSED_R1.fastq.gz'
R2=$sample'-QUALITY_PASSED_R1.fastq.gz'
flash $R1 $R2 -x 0.05 -m 20 -M 150 -z -o '$sample' 2>&1 | tee $sample'-flash.log'

sample="TARA_210"
R1=$sample'-QUALITY_PASSED_R1.fastq.gz'
R2=$sample'-QUALITY_PASSED_R1.fastq.gz'
flash $R1 $R2 -x 0.05 -m 20 -M 150 -z -o '$sample' 2>&1 | tee $sample'-flash.log'

# # ----- download MIME metagenomic datasets -----

# /users/home/cat3/projects/mime-mg/data/tara
mkdir -p $WD/data/mime
cd $WD/data/mime

echo -e "ERR5001722\nERR5005336\nERR5005972" > acc.txt

cat $WD/data/mime/acc.txt | while read -r acc
do
echo $acc
l=$(expr length $acc)

if [[ "$l" == 10 ]]; then
wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/00'${acc:9:10}'/'$acc'/'$acc'_1.fastq.gz'
wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/00'${acc:9:10}'/'$acc'/'$acc'_2.fastq.gz'
else
wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/'$acc'/'$acc'_1.fastq.gz'
wget 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'${acc:0:6}'/'$acc'/'$acc'_2.fastq.gz'
fi
done

echo -e "sample\tr1\tr2
MIME_001\tERR5001722_1.fastq.gz\tERR5001722_2.fastq.gz
MIME_002\tERR5005336_1.fastq.gz\tERR5005336_2.fastq.gz
MIME_003\tERR5005972_1.fastq.gz\tERR5005972_2.fastq.gz" > samples.txt

iu-gen-configs samples.txt -o $WD/data/mime
iu-filter-quality-minoche $WD/data/mime/MIME_001.ini --ignore-deflines
iu-filter-quality-minoche $WD/data/mime/MIME_002.ini --ignore-deflines
iu-filter-quality-minoche $WD/data/mime/MIME_003.ini --ignore-deflines

gzip *.fastq

cd $WD/data/mime
sample="MIME_001"
R1=$sample'-QUALITY_PASSED_R1.fastq.gz'
R2=$sample'-QUALITY_PASSED_R1.fastq.gz'
flash $R1 $R2 -x 0.05 -m 20 -M 150 -z -o $sample 2>&1 | tee $sample'-flash.log'


# ---- start with Synechococcus -----

mkdir -p $WD/data/METABAT__194-Synechococcus
cd $WD/data/METABAT__194-Synechococcus

wget 'ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737978/METABAT__194-surface-contigs.fa.gz'
gunzip METABAT__194-surface-contigs.fa.gz

anvi-gen-contigs-database -f METABAT__194-surface-contigs.fa -o METABAT__194-contigs.db -n Synechococcus -T 10

# Input FASTA file .............................: /users/work/cat3/projects/mime-mg/data/METABAT__194-Synechococcus/METABAT__194-surface-contigs.fa
# Name .........................................: Synechococcus
# Description ..................................: No description is given
# Num threads for gene calling .................: 10
#
# Finding ORFs in contigs
# ===============================================
# Genes ........................................: /tmp/tmpls98m1al/contigs.genes
# Amino acid sequences .........................: /tmp/tmpls98m1al/contigs.amino_acid_sequences
# Log file .....................................: /tmp/tmpls98m1al/00_log.txt
#
# CITATION
# ===============================================
# Anvi'o will use 'prodigal' by Hyatt et al (doi:10.1186/1471-2105-11-119) to
# identify open reading frames in your data. When you publish your findings,
# please do not forget to properly credit their work.
#
# Result .......................................: Prodigal (v2.6.3) has identified
#                                                 1809 genes.
# CONTIGS DB CREATE REPORT
# ===============================================
# Split Length .................................: 20,000
# K-mer size ...................................: 4
# Skip gene calling? ...........................: False
# External gene calls provided? ................: False
# Ignoring internal stop codons? ...............: False
# Splitting pays attention to gene calls? ......: True
# Contigs with at least one gene call ..........: 187 of 187 (100.0%)
# Contigs database .............................: A new database,
#                                                 METABAT__194-contigs.db, has
#                                                 been created.
# Number of contigs ............................: 187
# Number of splits .............................: 189
# Total number of nucleotides ..................: 1,553,062
# Gene calling step skipped ....................: False
# Splits broke genes (non-mindful mode) ........: False
# Desired split length (what the user wanted) ..: 20,000
# Average split length (what anvi'o gave back) .: 22,249

anvi-run-kegg-kofams -c METABAT__194-contigs.db -T 10

anvi-estimate-metabolism -c METABAT__194-contigs.db --only-complete

sample=TARA_158
TARADIR=/users/home/cat3/projects/mime-mg/data/tara
TARASMP=/users/home/cat3/projects/mime-mg/data/tara/$sample.extendedFrags.fastq.gz
echo $TARASMP
coverm genome --single $TARASMP --reference METABAT__194-surface-contigs.fa --single-genome -o $sample'-rel-abund.tsv'

# ---- check coverage (recruitement) of reads by the mags -----

mkdir -p /users/home/cat3/projects/mime-mg/data/mags
cd /users/home/cat3/projects/mime-mg/data/mags

# download all the mags
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737948/METABAT__121-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737976/METABAT__181-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737953/METABAT__133-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738037/METABAT__9-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737965/METABAT__157-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738014/METABAT__278-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737973/METABAT__178-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737975/METABAT__180-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737991/METABAT__22-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738010/METABAT__274-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738027/METABAT__67-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737956/METABAT__136-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737943/METABAT__116-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737978/METABAT__194-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737993/METABAT__224-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1737959/METABAT__142-surface-contigs.fa.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738030/METABAT__71-surface-contigs.fa.gz
gunzip *.fa.gz

for f in *.fa; do
    mv -- "$f" "${f%.fa}.fna"
done

sample=TARA_158
TARADIR=/users/home/cat3/projects/mime-mg/data/tara
TARASMP=/users/home/cat3/projects/mime-mg/data/tara/$sample.extendedFrags.fastq.gz
echo $TARASMP

coverm genome --genome-fasta-directory /users/home/cat3/projects/mime-mg/data/mags/ --dereplicate
    --single $TARASMP --single-genome -o $sample'-rel-abund.tsv'

sample=MIME_002
MIMEDIR=/users/home/cat3/projects/mime-mg/data/mime
MIMESMP='/users/home/cat3/projects/mime-mg/data/mime/'$sample'.extendedFrags.fastq.gz'
echo $MIMESMP

coverm genome --genome-fasta-directory /users/home/cat3/projects/mime-mg/data/mags/ --dereplicate --single $MIMESMP -o $sample'-rel-abund.tsv'
