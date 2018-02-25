#!/bin/bash

# https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1/$FILE_ID/$SUBDIR/$OUTPUT_FILE

FILE_NAME=aseg.stats

while IFS= read -r FILE_ID
do
    wget https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/freesurfer/5.1/$FILE_ID/stats/$FILE_NAME -O out/$FILE_ID.$FILE_NAME
done < FILE_ID.txt
