# Tips

### Units

Due to chemistry database PelePhysics, Cerisse uses cgs system.  
(is a pain, minor price for a bit faster chemistry and integration)

CHECKS

### Install yt

To load data into python for plotting, etc

```
$ python -m pip install --user yt
```
or using conda environments

```
$ conda install --channel conda-forge yt
```

Beware of the combination  Python **3.9** and **yt**. It is recommended to use newer or older Python, **>3.11** recommended.
Check [yt Website](https://yt-project.org) for details and tutorials.

### Use cerisse script

The **cerisse** script can help to clean directories, create backups, create visit files for movies quickly. It is useful for day to day running and avoid file creep-up.
Is installed in `bin`
and can be used by  

```
$ /bin/cerisse [options]
```

or by setting the path in `.bashrc` or `.profile` by adding the line at the end

```
PATH="/home/snm/work/cerisse/bin:$PATH"
```

where the `/home/snm/work/cerisse`corresponds to the installation dir


### Running on Imperial's HPC

To run on Imperial's CX2/3, you need to load the following modules

```bash
$ module load tools/prod iimpi/2020a
```

A typical job script looks like this:

```bash
#PBS -l select=2:ncpus=128:mem=256gb
#PBS -l walltime=24:00:00
#PBS -N name_of_the_job

module load tools/prod iimpi/2022a

export FI_MLX_IFACE=eth0
export I_MPI_HYDRA_BOOTSTRAP="rsh"
export I_MPI_HYDRA_BOOTSTRAP_EXEC="/opt/pbs/bin/pbs_tmrsh"
export I_MPI_HYDRA_BRANCH_COUNT=0

cd $PBS_O_WORKDIR

mpirun ./Cerisse3d.gnu.MPI.ex inputs
```

NOTE: There is an issue with MPI version >2020b on CX3 that causes deadlocks when writing chk or plt files.

#### Legacy (intel-2019.8.254)

Below is for running with Intel MPI 2019, which is unrecommended by RCS but still works.

```bash
$ module load gcc/11.2.0 mpi/intel-2019.8.254
```

If you see an error relating to `icpc: command not found`, that is because the `mpicxx` wrapper uses Intel's compiler by default. To change this to `g++`, you need to add the following command before running `make`, as

```bash
$ I_MPI_CXX=g++ make
```

