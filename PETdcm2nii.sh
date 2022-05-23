#!/usr/bin/env bash 

# mkdir -p ./demo1nii
# mkdir -p ./demo2nii

# dcm2niix -o ./demo1nii /Users/mac/Downloads/Epilepsy_group_3DT1WI/test1



for subj in `ls ./data/`
do
	mkdir -p ./data/${subj}/PET
	dcm2niix -o ./data/${subj}/PET  /home/wane/Documents/MRI阴性PET图像/data/${subj}/Static_Brain_3D_MAC/
	rm -rf ./data/${subj}/PET/*.json
done

