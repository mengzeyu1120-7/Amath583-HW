#!/bin/bash
#SBATCH --job-name=image_filter
#SBATCH --output=timing_results.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
#SBATCH --partition=gpu-rtx6k
#SBATCH --account=amath

echo "Compiling program..."
source cmpl.png-gs

echo "Threads,SeqTime(ms),ThreadedTime(ms)" > timing_data.csv
echo "Validation Results" > validation_results.txt

for threads in 1 2 3 4; do
    echo -e "\n\n===== Testing with $threads threads =====" | tee -a timing_results.txt
    
    ./xcpp-png-gs input_data.png output_data.png $threads | tee -a timing_results.txt
    
    seq_time=$(grep "Sequential" timing_results.txt | tail -1 | awk '{print $4}')
    thread_time=$(grep "Threaded" timing_results.txt | tail -1 | awk '{print $6}')
    
    echo "$threads,$seq_time,$thread_time" >> timing_data.csv
    
    echo -e "\nValidation for $threads threads:" | tee -a validation_results.txt
    ./xcompare_png output_data_seq.png output_data_thread.png | tee -a validation_results.txt
    
    if [ $threads -eq 4 ]; then
        cp output_data_thread.png threaded_4_result.png
        echo "Saved 4-thread result to threaded_4_result.png" | tee -a timing_results.txt
    fi
done

echo -e "\n\n===== Final Results =====" | tee -a timing_results.txt
echo "Timing data saved to timing_data.csv" | tee -a timing_results.txt
echo "Validation results saved to validation_results.txt" | tee -a timing_results.txt
echo "4-thread image saved to threaded_4_result.png" | tee -a timing_results.txt
cat timing_data.csv | tee -a timing_results.txt