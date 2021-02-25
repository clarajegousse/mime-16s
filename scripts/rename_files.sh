RUN_NUM=20190503_0073

cd /users/home/cat3/projects/mime-16s/data/miseq/$RUN_NUM

echo "Sample, R1, R2" > samples.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)

  rev_filename=$(echo $i | sed 's/R1/R2/')

  sample_name=$(echo $i | cut -f 1 -d "_")

  fwd_new_filename=$(echo $sample_name"_R1.fastq.gz" )
  rev_new_filename=$(echo $sample_name"_R2.fastq.gz" )

  cp $fwd_filename $fwd_new_filename
  cp $rev_filename $rev_new_filename

  echo "$sample_name, $fwd_new_filename, $rev_new_filename" >> samples.csv
done

cd /users/home/cat3/projects/mime-16s
mkdir -p /users/home/cat3/projects/mime-16s/results/$RUN_NUM/cutadapt
