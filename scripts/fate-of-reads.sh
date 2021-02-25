#!/bin/bash

MISEQ_RUN=20190508_0074

# raw reads from the miseq
RAW_READS_DIR=/users/home/cat3/projects/mime-16s/data/miseq/$MISEQ_RUN

# reads after adapter removal
CLEAN_READS_DIR=/users/home/cat3/projects/mime-16s/results/$MISEQ_RUN/cutadapt

# after dada2 trimming
TRIMMED_READS_DIR=/users/home/cat3/projects/mime-16s/results/$MISEQ_RUN/dada2-filter

cd $RAW_READS_DIR; touch read_counts.csv
echo "Sample, R1, R2" > read_counts.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filename=$(echo $i | sed 's/R1/R2/')
  sample_name=$(echo $i | cut -f 1 -d "_")

  fwd_num_raw=$(zcat $i | echo $((`wc -l`/4)))
  rev_num_raw=$(zcat $rev_filename | echo $((`wc -l`/4)))

  echo "$sample_name, $fwd_num_raw, $rev_num_raw" >> read_counts.csv
done

cd $CLEAN_READS_DIR; touch read_counts.csv
echo "Sample, R1, R2" > read_counts.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filename=$(echo $i | sed 's/R1/R2/')
  sample_name=$(echo $i | cut -f 1 -d "_")

  fwd_num_raw=$(zcat $i | echo $((`wc -l`/4)))
  rev_num_raw=$(zcat $rev_filename | echo $((`wc -l`/4)))

  echo "$sample_name, $fwd_num_raw, $rev_num_raw" >> read_counts.csv
done

cd $TRIMMED_READS_DIR; touch read_counts.csv
echo "Sample, R1, R2" > read_counts.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filename=$(echo $i | sed 's/R1/R2/')
  sample_name=$(echo $i | cut -f 1 -d "_")

  fwd_num_raw=$(zcat $i | echo $((`wc -l`/4)))
  rev_num_raw=$(zcat $rev_filename | echo $((`wc -l`/4)))

  echo "$sample_name, $fwd_num_raw, $rev_num_raw" >> read_counts.csv
done
