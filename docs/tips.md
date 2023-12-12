# Tips

### Units

Due to chemistry database PelePhysics, Cerisse uses cgs system.  
(is a pain, minor price for a bit faster chemistry and integration)

### Install yt

To load data into python for plotting, etc

```
$ python -m pip install --user yt
```
or using conda environments

```
$ conda install --channel conda-forge yt
```

Beware of Python **3.9** and **yt** use newer or older, **>3.11** recommended.
Check [yt Website](https://yt-project.org) for details and tutorials.

### Use cerisse_help script

Script file to clean directories, create backups, create visit file for movies.
Useful for day to day running and avoid file creep-up.
Is installed in `EB_CNS/Tools/scripts`
and can be used by  

```
$ ../../Tools/scripts/cerisse_help [options]
```

or by setting the path in `.bashrc` or `.profile` by adding the line at the end

```
PATH="/home/snm/work/cerisse/EB_CNS/Tools/scripts:$PATH"
```

where the `/home/snm/...`coresponds to the installation dir


### Defaults

What are numeric defaults?


