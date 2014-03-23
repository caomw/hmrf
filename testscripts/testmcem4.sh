myhost=`hostname`
session_dir=$1
initlabelmap=$2

test_dir=testmcem4
if [ ! -d ${session_dir}/${test_dir} ]
then 
    mkdir ${session_dir}/${test_dir}
fi

~/projects/mcem_lemon/build_${myhost}/mcem_lemon --burnin 20 --numSamples 50 --emiter 30 --beta 10 --gamma 1 -k 7 --inittemp 3 --finaltemp 3  -i ${initlabelmap}  -f  ${session_dir}/func_sphere.nii.gz --cumusample ${session_dir}/${test_dir}/cumusample.nii.gz --rsample ${session_dir}/${test_dir}/rsample.nii.gz -v 3

~/projects/group_mrf_lemon/build_${myhost}/sample2label -i $session_dir/$test_dir/cumusample.nii.gz -o $session_dir/$test_dir/labelmap.nii.gz
~/projects/group_mrf_lemon/build_${myhost}/copyheader -s $initlabelmap -i $session_dir/$test_dir/labelmap.nii.gz -o $session_dir/$test_dir/labelmap.nii.gz
