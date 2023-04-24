bsub -n 33 -W 24:00 -J train bash -c 'mpirun ./app_main --nEnvironments 4 --nMasters 1 --nThreads 1 --workerProcessesPerEnv 8'
