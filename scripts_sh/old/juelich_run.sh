#!/bin/bash
#SBATCH --job-name=slim_1
#SBATCH --nodes=1
#SBATCH --ntasks=46
#SBATCH --account=eiri_e_unipa2       # Your budgeting account
#SBATCH --partition=g100_usr_prod   # Set the partition to run on
#SBATCH -o mouse_%j.out            # File to which STDOUT will be written 
#SBATCH -e mouse_%j.err            # File to which STDERR will be written 
#SBATCH --time=04:00:00

export HDF5_USE_FILE_LOCKING=FALSE

source $WORK/SW/Environments/load_modules.sh
#source $WORK/SW/Environments/cereb_env/bin/activate 
source $WORK/SW/nest-3.7/bin/nest_vars.sh
source $HOME/bsb-env/bin/activate
#mpirun -n 46 bsb compile --clear -v 3 configurations/mouse/nest/basal_vitro.yaml
#mpirun -n 48 python run_reconstruct.py
#mpirun -n 46 python run_sim.py mouse_cerebellum.hdf5

# here
source $XXX/sim_setup.sh