#!/bin/bash

## with this bash script, you can generate a bash file wich can perform
## bach processing for test images
## Note: you may need to slightly modify this script

data_dir='<your dir>' # please input your private dir
imgs=`find ${data_dir} -name "*_A.png"`
for img in ${imgs}; do
    echo ${img}
    filename=`basename ${img}`
    echo "python infer.py -c jshdr -i ./data/test/${filename}" >> test_list.sh
done
