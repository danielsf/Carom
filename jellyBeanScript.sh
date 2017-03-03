make jellyBean_maps
#make analysis

out_dir=output/workspace/

# seeds 99 66 125 6475

#for ss in 626 694 762 1068 6475 66 125
#for ss in 66 694 6475 762 1068 125 626
#for ss in 66 694 1068 626 6475 762
for ss in 11 452 831 642 53 19 90 52 112 7742 5345 1082
do
    out_name=${out_dir}jellyBean_d12_s${ss}_output.sav
    
    ./bin/jellyBean_maps -d 12 -c 21.03 -p 0.95 -n 400000 -s ${ss} \
    -t ${out_dir}jellyBean_d12_s${ss}_timing.sav \
    -o ${out_name} -i 24 >& ${out_dir}jellyBean_s${ss}_stdout.txt
    
    #./bin/analysis -i ${out_name} -o output/scratch/jellyBean_d12_s${ss}_processed \
    #-d 12 -c 21.0 -x 0 1 -x 6 9
done


#./bin/jellyBean_maps -d 12 -c 21.0 -p 0.95 -n 300000 -s 99 -t output/scratch/jellyBean_d12_s99_timing.sav -o output/scratch/jellyBean_d12_s99_output.sav

#./bin/jellyBean_maps -d 12 -c 21.0 -p 0.95 -n 300000 -s 66 -t output/scratch/jellyBean_d12_s66_timing.sav -o output/scratch/jellyBean_d12_s66_output.sav

#./bin/jellyBean_maps -d 12 -c 21.0 -p 0.95 -n 300000 -s 125 -t output/scratch/jellyBean_d12_s125_timing.sav -o output/scratch/jellyBean_d12_s125_output.sav

#./bin/jellyBean_maps -d 12 -c 21.0 -p 0.95 -n 300000 -s 6475 -t output/scratch/jellyBean_d12_s6475_timing.sav -o output/scratch/jellyBean_d12_s6475_output.sav
