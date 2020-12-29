# Generate impact plots 2020/12

## generate SoLID pseudo data
* compile code

make O=analysis_neutron

make O=analysis_proton

* get pseudodata

./analysis_neutron 1 //generate bins

./analysis_neutron 2 //generate projection root files

./analysis_neutron 3 //convert root files into text table

similar for the proton case

## fit

* prepare pseudodata sets

copy generated SoLID pseudodata into ./data and ./datacollins

./prepare sivers //prepare pseudodata sets for Sivers fit

./prepare collins //prepare pseudodata sets for Collins fit

* run the fit

./fitsivers.py [opt] //see opt values at the end part of the code

./fitcollins.py [opt]

* plot

use jupyter-notebook open the plot notebook to make plots


