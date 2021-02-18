cd /users/home/cat3/projects/mime-16s/data/miseq/20201112_0114

echo "Sample, R1, R2" > samples.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filnemae=$(echo $i | sed 's/R1/R2/')

  samplename=$(echo $i | cut -f 1 -d "_")

  mv $fwd_filename $fwd_new_filename
  mv $rev_filename $rev_new_filename

  fwd_new_filename=$(echo $samplename"_R1.fastq.gz" )
  rev_new_filename=$(echo $samplename"_R2.fastq.gz" )

  echo "$samplename, $fwd_new_filename, $rev_new_filename" >> samples.csv
done

cd /users/home/cat3/projects/mime-16s
mkdir -p /users/home/cat3/projects/mime-16s/20200416_0101/cutadapt
