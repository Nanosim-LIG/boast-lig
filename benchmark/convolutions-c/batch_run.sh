#!/bin/bash

# See "man qsub" for an extensive description of options
# See http://www.clusterresources.com/torquedocs21/2.1jobsubmission.shtml for a description of resource requests

# Output files will be named "job-name.{o,e}[0-9]*"
#PBS -N convolution-test
# Get two full nodes, one process per node
# Please use "procs=N" if you do not care about exclusive node acces
#PBS -l nodes=1
# Queue name (see "qstat -q" for a list)
#PBS -q batch
# Expected execution time (optional; helps the scheduler improve cluster utilization)
#PBS -l walltime=00:25:00
# Do not re-run the job if it fails
#PBS -r n
# Set the working directory
#PBS -d /home/jcronsioe/bigdft-1.6.0/src/convolutions-c 
#PBS -k eo
 
# Be nice to the system if we try to allocate more physical memory than it is
# available (otherwise the system could partially crash).
echo 1000 >> /proc/self/oom_score_adj

# Make the given program span all the allocated nodes for the job.
#
# We also need to make mpiexec children behave nicely (w.r.t. OOM).
mpiexec sh -c 'echo 1000 >> /proc/self/oom_score_adj; echo "I am process $$ (son of $PPID) at `hostname` oom=`cat /proc/self/oom_score`"'

sh run.sh 
