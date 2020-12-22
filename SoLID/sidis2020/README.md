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
