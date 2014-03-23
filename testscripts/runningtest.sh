myhost=`hostname`
session_dir=$1
initgrpimage=$2
test_dir=$3

if [ ! -d ${session_dir}/${test_dir} ]
then
    mkdir ${session_dir}/${test_dir}
    mkdir ${session_dir}/${test_dir}/cumusamples
    mkdir ${session_dir}/${test_dir}/rsamples
    mkdir ${session_dir}/${test_dir}/labelmap
fi

~/projects/group_mrf_lemon/build_${myhost}/groupmrf --nthreads 300 --sweepperthread 1 --burnin 10 --numSamples 10 --emiter 2 --alpha 1 --beta 0.2 --gamma 1 -k 7 --estprior --inittemp 1.0 --finaltemp 1.0 --initsame -i ${initgrpimage}  -f  ${session_dir}/allfunc_normalized --cumusampledir ${session_dir}/${test_dir}/cumusamples --rsampledir ${session_dir}/${test_dir}/rsamples -v 2

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
