# Calculate energy grids from density grids 
'Manual' approach to derivate solvent-probes binding free energies from density grids generated in MDMix.

This functionality is usually covered normally in mdmix, but due to the rigidity of mdmix projects it might be useful.

The algorithm I used varies from MDMix. pyMDMix version of the energy calculation can be found [here](http://mdmix.sourceforge.net/analysis-guide/#energy).
Basically it obtains the estimated density from the concentration and corrects the energy by the differences in volume from out simulation and the standart concentrations (1M). 

The algorithm I propose here is independent of the concentration and number of frames stored, as it calculates the global appearence of probes and estimates the uniform baseline. Then it estimates the binding free energy from the difference with this density and the simulated one. I have not impremented the volume correction.

##1. Usage

To obtain the energy grid from the three replicas use _dgrids2egrid.py_
```{bash}
python lib/dgrids2egrid.py dgrid1.dx dgrid2.dx dgrid3.dx output_egrid_AVG.dx
```

To obtain the energy from a single density grid use
```{bash}
python lib/dgrid2egrid.py dgrid.dx output_egrid.dx
``` 
