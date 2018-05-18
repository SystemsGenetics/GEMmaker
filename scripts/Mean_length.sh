#!/bin/bash
#This script calculates average length of fastq files.

#This brings in the parameters of the name of the file, and the minlen percentage dictated by the user
srx=$1
minlen=$2
#cd "../${srx}"
total=0

#This if statement checks if the data is single or paired data, and checks length accordingly
#This script returns 1 number, which can be used for the minlen in trimmomatic
if [ -e ${srx}_1.fastq ] && [ -e ${srx}_2.fastq ]; then
  for fastq in "${srx}_1.fastq" "${srx}_2.fastq"; do
    a=`awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' $fastq \
    | sort \
    | awk '{ print $0, $1*$2}' \
    | awk '{ SUM += $3 } { SUM2 += $2 } END { printf("%.0f\n", SUM / SUM2 * '$minlen')} '`
  #echo $minlen
  total=$(( ($a + $total) ))
  done
  total=$(( $total / 2 ))
  echo $total

elif [ -e ${srx}_1.fastq ]; then
  awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ${srx}_1.fastq \
    | sort \
    | awk '{ print $0, $1*$2}' \
    | awk '{ SUM += $3 } { SUM2 += $2 } END { printf("%.0f\n", SUM / SUM2 * '$minlen')} '
fi

