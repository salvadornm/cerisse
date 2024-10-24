# Chemistry



## PelePhysics


Uses the [PelePhysics](https://pelephysics.readthedocs.io/en/latest/index.html)

In the `GNUmakefile` file you need to need to set-up the options

```bash 
# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := BurkeDryer
```

and make sure the following lines exists in that file

```
USE_PELEPHYSICS = TRUE
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)
```


### Available Mechanisms


### Available EOS Models

Using PelePhysics there are three EOS models: 

- A simple GammaLaw model for a single component perfect gas

- An ideal gas mixture model labeled **Fuego**

- The Soave-Redlich-Kwong cubic equation of state for a general mixture of species 

### Generate a New mechanism

For all the available mechanisms, a Cantera yaml format is provided.
If CHEMKIN fiels are present
Pelephysics rely on Cantera’s ***ck2yaml*** utility to convert CHEMKIN files to the Cantera yaml format

Check [PelePhysics Tutorial ](https://pelephysics.readthedocs.io/en/latest/EOS.htmlls\)


