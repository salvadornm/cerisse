---
icon: info
---

# Tips

Small tips on installation of auxiliary files and scripts

## yt

[yt](https://yt-project.org) is an open-source python package for analyzing and visualizing volumetric data. It us highly recomended if Python is used for postprocessing and analysis of simulations data. To install

```bash
$ python -m pip install --user yt
```

or using conda environments

```bash
$ conda install --channel conda-forge yt
```

Beware of the combination Python **3.9** and **yt**. It is recommended to use newer Python, **>3.11** recommended. Check [yt Website](https://yt-project.org) for details and tutorials. Most of the examples in the manual use it for quick analysis and plotting in a `plot.py` script.

In recent systems (for example Ubuntu 24), the use of virtual environments is enforced and the above lines will not work unless a virtual environment is used (or altenative ways created). For example, if the ```virtual``` environment exist, to activate it type:

```bash
$ source ~/virtual/bin/activate
```

You will need to do this to install python package, such as yt and markdown, required for postprocessing or managing the documentation.

## Cerisse script

The **cerisse** script can help to clean directories, create backups, create visit files for movies quickly. It is useful for day to day running and avoid file creep-up. The script is installed in `bin` and can be used by

```bash
$ ../bin/cerisse [options]
```

or by setting the path in `.bashrc` or `.profile` by adding the line at the end

```
PATH="/home/snm/work/cerisse/bin:$PATH"
```

where `/home/snm/work/cerisse` corresponds to the installation directory. Current options are **clean** (_all_) to remove old results and temporary directories. **visit** to prepare data for input into Visit, **plot**, which will crete a snaphsot of the density, **backup** (_dirname_) will save simualytion data and main files into a backup directory

## Documentation Editing

Althoigh the main documentation is onine in Gitbooks. It is possible to generate the documentation by locally typing in the parent directory: `$ mkdocs serve` and the point the browser to [127.0.0.1.8000](http://127.0.0.1:8000) in your browser. The formatting of this documentation may be off (inline equations in particular) and not all images wikk appear.

You may need to install the `python-markdown-math` extension for rendering equations and the `markdown-callouts` extension for correctly displaying the warning and note blocks. For help editing the documentation visit [mkdocs.org](https://www.mkdocs.org).

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
