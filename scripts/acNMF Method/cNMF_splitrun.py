import csv
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scanpy as sc
from IPython.display import Image

np.random.seed(14)

########################
#Set cNMF run parameters
########################

numiter = 200
numhvgenes = 2000

K = ' '.join([str(i) for i in range(2,61)] + [str(i) for i in range(70,210,10)])
K_int = [int(i) for i in K.split()]

numworkers = numiter*len(K_int)

output_directory = "/<path>/<to>/<output_dir>/"
run_name = 'cNMF_split0'

countfn = "cnmf_input_benchmark_split0.h5ad"

seed = 14

###############################
#Set up directories for LSF i/o
###############################

logdir = os.getcwd() + '/split0_logdir'
os.system('mkdir ' + logdir)
factorize_dir = logdir + '/factorize'
os.system('mkdir ' + factorize_dir)
combine_dir = logdir + '/combine'
os.system('mkdir ' + combine_dir)
consensus_dir = logdir + '/consensus'
os.system('mkdir ' + consensus_dir)

##############
#Preprocessing
##############

prepare_cmd = 'python /home/rchapple/software/cNMF/cnmf_v2.0.py prepare --output-dir %s --name %s -c %s -k %s --n-iter %d --total-workers %d --seed %d --numgenes %d --beta-loss frobenius' % (output_directory, run_name, countfn, K, numiter, numworkers, seed, numhvgenes)
print('Prepare command assuming parallelization with %d tasks:\n%s' % (numworkers, prepare_cmd))
os.system(prepare_cmd)

##############
#Factorization
##############

worker_index = ' '.join([str(x) for x in range(numworkers)])

#Set up the job submission array
#Each job runs all iterations of NMF for each rank in one instance of opening the cNMF_v2.0.py script
#The number of jobs that get submitted to HPC are equivalent to the number of ranks

start = [int(i) for i in range(0, numworkers-1, numiter)]
end = [int(i) for i in range(numiter-1, numworkers, numiter)]
nmf_job_data = {'K':pd.Series(K_int), 'Start':pd.Series(start), 'End':pd.Series(end)}
nmf_job_submission = pd.DataFrame(nmf_job_data)

for x in nmf_job_submission.index:
    factorize_cmd = "bsub -P RC -J split0_factorize -R \"span[hosts=1] rusage[mem=10GB]\" -oo %s -eo %s 'python /home/rchapple/software/cNMF/cnmf_v2.0.py factorize --output-dir %s --name %s --jobstart %d --jobend %d'" % (factorize_dir, factorize_dir, output_directory, run_name, nmf_job_submission['Start'][x], nmf_job_submission['End'][x])
    print('Factorize command to run factorizations for rank = %d across all iterations' % (nmf_job_submission['K'][x]))
    os.system(factorize_cmd)

wait = "bwait -w 'ended(split0_factorize)'"
os.system(wait)

#################################
#Combine factorization replicates
#################################

combine_cmd = "bsub -P RC -J split0_combine -R\"rusage[mem=400GB]\" -q large_mem -oo %s -eo %s 'python /home/rchapple/software/cNMF/cnmf_v2.0.py combine --output-dir %s --name %s'" % (combine_dir, combine_dir, output_directory, run_name)
print(combine_cmd)
print('\n')
os.system(combine_cmd)

wait_combine = "bwait -w 'ended(split0_combine)'"
os.system(wait_combine)


##########################################################
#Cluster and plot consensus programs and KNN outlier plots
##########################################################

from itertools import chain

density_threshold = 2.00

for x in chain(range(2,61), range(70,210,10)):
    selected_K = x
    print("Selected_K =", selected_K, "\n")
    consensus_cmd = "bsub -P RC -J split0_consensus -R\"rusage[mem=400GB]\" -q large_mem -oo %s -eo %s 'python /home/rchapple/software/cNMF/cnmf_v2.0.py consensus --output-dir %s --name %s --local-density-threshold %.2f --components %d --show-clustering'" % (consensus_dir, consensus_dir, output_directory, run_name, density_threshold, selected_K)
    print('Consensus command for K=%d:\n%s' % (selected_K, consensus_cmd))
    os.system(consensus_cmd)

density_threshold_str = ('%.2f' % density_threshold).replace('.', '_')

#Outlier filtering

density_threshold = 0.10                                                        

for x in chain(range(2,61), range(70,210,10)):                                                       
    selected_K = x                                                              
    print("Selected_K =", selected_K, "\n")                                     
    consensus_cmd = "bsub -P RC -J split0_consensus_outlier -R\"rusage[mem=400GB]\" -q large_mem -oo %s -eo %s 'python /home/rchapple/software/cNMF/cnmf_v2.0.py consensus --output-dir %s --name %s --local-density-threshold %.2f --components %d --show-clustering'" % (consensus_dir, consensus_dir, output_directory, run_name, density_threshold, selected_K)
    print('Consensus command for K=%d:\n%s' % (selected_K, consensus_cmd)) 
    os.system(consensus_cmd)
                                                        
density_threshold_str = ('%.2f' % density_threshold).replace('.', '_')
