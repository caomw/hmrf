myhost=`hostname`
# same with test2b but beta is weighted by spatial distance.
session_dir=$1
initgrpimage=$2
test_dir=test2f

if [ ! -d ${session_dir}/${test_dir} ]
then
    mkdir ${session_dir}/${test_dir}
    mkdir ${session_dir}/${test_dir}/cumusamples
    mkdir ${session_dir}/${test_dir}/rsamples
    mkdir ${session_dir}/${test_dir}/labelmap
fi

~/projects/group_mrf_lemon/build_${myhost}/groupmrf --nthreads 300 --sweepperthread 1 --burnin 20 --numSamples 100 --emiter 30 --alpha 1 --beta 3 --weightbetadist --gamma 1 -k 7 --inittemp 0.1 --finaltemp 0.1 --initsame -i ${initgrpimage}  -f  ${session_dir}/allfunc_normalized --cumusampledir ${session_dir}/${test_dir}/cumusamples --rsampledir ${session_dir}/${test_dir}/rsamples

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
