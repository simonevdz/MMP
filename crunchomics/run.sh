#!/bin/bash
n=$1

source /zfs/omics/software/v_envs/madpy3/bin/activate

for i in $(seq 1 $n)
do
  sbatch RNA-seq.sh $i $n
done



