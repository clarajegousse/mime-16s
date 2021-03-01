MIME_DIR=/users/home/cat3/projects/mime-16s
RUN_NUM=20190503_0073

cd $MIME_DIR
cp -r ../mime/data/miseq/$RUN_NUM data/miseq/

cd $MIME_DIR/data/miseq/$RUN_NUM

echo "Sample, R1, R2" > samples.csv
for i in `ls *R1*.fastq.gz`; do

  fwd_filename=$(echo $i)
  rev_filename=$(echo $i | sed 's/R1/R2/')
  sample_name=$(echo $i | cut -f 1 -d "_" | sed 's/-Archaea/-ARK/')
  #echo $sample_name

  fwd_new_filename=$(echo $sample_name"_R1.fastq.gz" )
  rev_new_filename=$(echo $sample_name"_R2.fastq.gz" )

  # cp $fwd_filename $fwd_new_filename
  # cp $rev_filename $rev_new_filename

  mv $fwd_filename $fwd_new_filename
  mv $rev_filename $rev_new_filename

  echo "$sample_name, $fwd_new_filename, $rev_new_filename" >> samples.csv
done

cd $MIME_DIR
mkdir -p /users/home/cat3/projects/mime-16s/results/$RUN_NUM/cutadapt

# git stash; git pull; chmod +x $MIME_DIR/scripts/*
