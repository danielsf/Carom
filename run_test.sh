make curved_12d

out_dir=$1

echo $out_dir

for ss in 771 11 831 452 642 53 19 52 112 90
do
    nohup nice ./bin/curved_12d -s ${ss} -n 600000 \
    -o ${out_dir}/output_s${ss}.txt \
    -t ${out_dir}/timing_s${ss}.txt >& ${out_dir}/stdout_${ss}.txt &

done
