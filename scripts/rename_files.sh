for i in `ls *R1*.fastq.gz`; do
  filename=$(echo $i)
  samplename=$(echo $i | cut -f 1 -d "_")
  newfilename=$(echo $samplename".1.fastq.gz" )
  mv $filename $newfilename
done

for i in `ls *R2*.fastq.gz`; do
  filename=$(echo $i)
  samplename=$(echo $i | cut -f 1 -d "_")
  newfilename=$(echo $samplename".2.fastq.gz" )
  mv $filename $newfilename
done
