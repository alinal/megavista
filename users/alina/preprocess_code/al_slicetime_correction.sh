#!/bin/bash
# al_slicetime_correction.sh
#

func_dir='/Volumes/Plata1/DorsalVentral/fmri/WC030311/WC030311_nifti'
newAdd='_stc'

## 6. Generate nifit file file (for use later)
echo "Running 3dTshift to get new time series"

cd $func_dir
output=$(ls *mcf.nii*)

for file in $output
do
    newFile=${file%%.*}
    newName=$newFile$newAdd
    3dTshift -verbose  -slice 16  -TR 2.0s  -tpattern seq-z  -prefix ${func_dir}/${newName}.nii.gz  ${func_dir}/${file}
done