#!/bin/bash

# aseg.stats
FILE_NAME=aseg.stats

mkdir $FILE_NAME

while IFS= read -r FILE_ID
do
    wget https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1/$FILE_ID/stats/$FILE_NAME -O $FILE_NAME/$FILE_ID.$FILE_NAME
done < FILE_ID.txt


# lh.aparc.stats
FILE_NAME=lh.aparc.stats

mkdir $FILE_NAME

while IFS= read -r FILE_ID
do
    wget https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1/$FILE_ID/stats/$FILE_NAME -O $FILE_NAME/$FILE_ID.$FILE_NAME
done < FILE_ID.txt


# rh.aparc.stats
FILE_NAME=rh.aparc.stats

mkdir $FILE_NAME

while IFS= read -r FILE_ID
do
    wget https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1/$FILE_ID/stats/$FILE_NAME -O $FILE_NAME/$FILE_ID.$FILE_NAME
done < FILE_ID.txt
