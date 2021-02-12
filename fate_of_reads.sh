

cd /users/home/cat3/projects/mime-16s/data/miseq/20190508_0074

echo "Sample, R1, R2" > read_counts.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filname=$(echo $i | sed 's/R1/R2/')
  samplename=$(echo $i | cut -f 1 -d "_")

  fwd_num_raw=$(zcat $i | echo $((`wc -l`/4)))
  rev_num_raw=$(zcat $rev_filname | echo $((`wc -l`/4)))

  echo "$samplename, $fwd_num_raw, $rev_num_raw" >> read_counts.csv

done
