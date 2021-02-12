cd /users/home/cat3/projects/mime-16s/data/miseq/20190508_0074

echo "Sample, R1, R2" > samples.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filnemae=$(echo $i | sed 's/R1/R2/')

  samplename=$(echo $i | cut -f 1 -d "_")

  mv $fwd_filename $fwd_new_filename
  mv $rev_filnemae $rev_new_filename

  fwd_new_filename=$(echo $samplename"_R1.fastq.gz" )
  rev_new_filename=$(echo $samplename"_R2.fastq.gz" )

  echo "$samplename, $fwd_new_filename, $rev_new_filename" >> samples.csv

done
