BATCH --partition=compute

#SBATCH --nodes=3

#SBATCH --ntasks-per-node=16

#SBATCH --time=5-00:00:00

#SBATCH --job-name=100nsT2POPC

#SBATCH --constraint=haswell


#SBATCH --mail-type=ALL

#SBATCH --mail-user=zoe.ford@pmb.ox.ac.uk  


module purge

module load gromacs/2018__single


. enable_arcus-b_mpi.sh

echo "MPI_HOSTS=${MPI_HOSTS}"


mpirun $MPI_HOSTS gmx_mpi mdrun -v -stepout 10000 -deffnm md 
