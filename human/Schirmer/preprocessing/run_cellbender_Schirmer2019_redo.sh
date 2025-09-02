#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/MS_Schirmer_2019

# get a list of all samples:
sample_dir="/dfs3b/swaruplab/shared_lab/cross-disorder/count_matrices/MS_Schirmer_2019/"
# samples=($(ls $sample_dir))


# samples that I need to re-do
# samples=("SRR9123032" "SRR9123033" "SRR9123038" "SRR9123039" "SRR9123040" "SRR9123042" "SRR9123043" "SRR9123046" "SRR9123048" "SRR9123050" "SRR9123051" "SRR9123052")
# expected_cells=(400 4000 400 400 400 400 400 400 400 400 400 400)

samples=("SRR9123033" "SRR9123038" "SRR9123039" )
expected_cells=(10000 10000 10000)

END=3
for i in $(seq 0 $END);
do
    echo "first redo"
    echo $i
    cur="${samples[$i]}"
    cells="${expected_cells[$i]}"

    # set input and output files:
    infile=$sample_dir$cur/counts_unfiltered/
    outfile=$sample_dir$cur/counts_unfiltered/cellbender.h5

    # run cellbender
    cellbender remove-background \
       --input $infile \
       --output $outfile \
       --expected-cells $cells \
       --total-droplets-included 25000 \
       --epochs 150 \
       --fpr 0.01 \
       --cuda

done





# samples that I need to re-do
samples=("SRR9123040" "SRR9123042" "SRR9123043" "SRR9123046" "SRR9123048" "SRR9123050" "SRR9123052")
expected_cells=(400 600 800 1000 1200 1400)

for cur in ${samples[@]};
do
  for cells in ${expected_cells[@]};
  do
      echo "$i $cells"

      mkdir $sample_dir$cur/counts_unfiltered/$cells

      # set input and output files:
      infile=$sample_dir$cur/counts_unfiltered/
      outfile=$sample_dir$cur/counts_unfiltered/$cells/cellbender.h5

      # run cellbender
      cellbender remove-background \
         --input $infile \
         --output $outfile \
         --expected-cells $cells \
         --total-droplets-included 15000 \
         --epochs 150 \
         --fpr 0.01 \
         --cuda

  done
done
