# create sample sheet for RNAseq data
# if it exists, ask user to remove it
if [ -f /projects/fs1/esfontjo/nhp/data/RNAseq/unzipped/sample_sheet.csv ]; then
  echo "sample_sheet.csv already exists. Please remove it and run again."

# if it doesn't create it
else
  printf "sample,fastq_1,fastq_2,strandedness\n" > /projects/fs1/esfontjo/nhp/data/RNAseq/unzipped/sample_sheet.csv
  # loop over all zip files in the directory
  find /projects/fs1/esfontjo/nhp/data/RNAseq -type f -name "*.zip" | while read line; do
    #  echo $line
    # remove path and extension
    sample=$(basename "$line" .zip)
    #  echo $sample
    # add line to file
    data_path="/projects/fs1/esfontjo/nhp/data/RNAseq/unzipped/"
    newline="${sample},${data_path}${sample}_1.fq.gz,${data_path}${sample}_2.fq.gz,unstranded"
    echo "$newline"  >> /projects/fs1/esfontjo/nhp/data/RNAseq/unzipped/sample_sheet.csv
  done
fi

