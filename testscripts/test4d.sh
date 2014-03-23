#!/bin/bash
# This test change beta, burnin and numClusters and save a slice in an png file.
BETA=$1
NUMCLUSTERS=$2
MAXSCANS=$3
SNAPINT=$4
pngdir=$5

myhost=`hostname`
session_dir=./preprocessed_nosmooth/session1_stage2b_less
initgrpimage=./Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz 
test_dir=test4d

rm $session_dir/$test_dir/rsamples/*
rm $session_dir/$test_dir/labelmap/*
rm $session_dir/$test_dir/slices/*
rm $session_dir/$test_dir/pngs/*

if [ ! -d $session_dir/$test_dir/$pngdir ]
then 
    mkdir $session_dir/$test_dir/$pngdir 
else
    rm $session_dir/$test_dir/$pngdir/*
fi 

~/projects/group_mrf_lemon/build_${myhost}/groupmrf --nthreads 100 --sweepperthread 1 --burnin $MAXSCANS --snapshotfreq $SNAPINT --numSamples 1 --emiter 1 --alpha 0 --beta $BETA --gamma 0 -k $NUMCLUSTERS --inittemp 1.0 --finaltemp 1.0 --initsame -i ${initgrpimage}  -f  ${session_dir}/allfunc_normalized --cumusampledir ${session_dir}/${test_dir}/cumusamples --rsampledir ${session_dir}/${test_dir}/rsamples

# post processing. 
cd $session_dir/$test_dir/rsamples
allfiles=`ls -d *`
cd -
for file in $allfiles
do
    ~/projects/group_mrf_lemon/build_${myhost}/sample2label -i $session_dir/$test_dir/rsamples/$file -o $session_dir/$test_dir/labelmap/$file
    fslroi $session_dir/$test_dir/labelmap/$file $session_dir/$test_dir/slices/$file 0 61 0 73 30 1
    bname=`basename $file .nii.gz`
    ~/packages/itkapp4.2_uv/ConvertBetweenFileFormats/ConvertBetweenFileFormats $session_dir/$test_dir/slices/$file $session_dir/$test_dir/pngs/${bname}.png unsigned_short
done

matlab -nodesktop -nosplash -nojvm -r "myind2cor('${session_dir}/${test_dir}/pngs', '${session_dir}/${test_dir}/${pngdir}', 8); exit;"