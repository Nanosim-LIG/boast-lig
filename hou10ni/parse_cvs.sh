#!/bin/bash

FILE=test.cvs
RES_FILE=res.dat

echo "=> Parsing file $FILE"
LINE_NUMBER=0
TIME_REF=0

printf "#Options\t Time\n" >> $RES_FILE
printf "\"\"\t 0.0000000\n" >> $RES_FILE

while read Line; do

	LINE_NUMBER="$(($LINE_NUMBER + 1))"
  OPTIONS=`echo $Line | cut -d',' -f1,2,3,4,5`
  TIME=$(echo $Line | cut -d',' -f6)
  printf -v TIME "%.7f" "$TIME"

	if [[ $LINE_NUMBER -eq 1 ]]; then
  	printf -v TIME_REF "%.7f" "$TIME"
  	printf "\"$OPTIONS\"\t $TIME\n " >> $RES_FILE
	fi

	if [[ $(echo "$TIME < $TIME_REF" | bc) -eq 1 ]]; then
  	printf "\"$OPTIONS\"\t $TIME\n " >> $RES_FILE
	fi


done < $FILE

echo "=> Done"
