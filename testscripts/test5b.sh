#!/bin/bash

# This test verify the chaging temperature.
myhost=`hostname`
session_dir=$1
initgrpimage=$2
test_dir=test5b

if [ ! -d ${session_dir}/${test_dir} ]
then
    mkdir ${session_dir}/${test_dir}
    mkdir ${session_dir}/${test_dir}/cumusamples
    mkdir ${session_dir}/${test_dir}/rsamples
    mkdir ${session_dir}/${test_dir}/labelmap
fi


~/projects/group_mrf_lemon/build_${myhost}/groupmrf --nthreads 300 --sweepperthread 1 --burnin 20 --numSamples 100 --emiter 50 --alpha 1 --beta 2 --gamma 1 -k 7 --inittemp 2 --finaltemp 0.1 --initsame -i ${initgrpimage}  -f  ${session_dir}/allfunc_normalized --cumusampledir ${session_dir}/${test_dir}/cumusamples --rsampledir ${session_dir}/${test_dir}/rsamples

# post processing. 

cd $session_dir/$test_dir/cumusamples
allfiles=`ls -d *`
cd -
labelmap="labelmap"
for file in $allfiles
do
    ~/projects/group_mrf_lemon/build_${myhost}/sample2label -i $session_dir/$test_dir/cumusamples/$file -o $session_dir/$test_dir/labelmap/$file
    ~/projects/group_mrf_lemon/build_${myhost}/copyheader -s $initgrpimage -i $session_dir/$test_dir/labelmap/$file -o $session_dir/$test_dir/labelmap/$file
done
