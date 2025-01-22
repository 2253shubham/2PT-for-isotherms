Quick steps to perform 2PT analysis

1. Generate velocity autocorrelation data from MD simulations (please refer "Computational details for methane adsorption in silicalite" section of our [publication](https://doi.org/10.1063/5.0099790), to understand how to set up the simulations and required parameters)	
2. Run [IAG_DoS_comp.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/HS_DoS_comp.py)/[HS_DoS_comp.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Hard_Sphere_approximation/scripts/HS_DoS_comp.py) on the velocity autocorrelation data obtained from MD simulations. 
i.	For IAG systems, you will need to run an extra code [refer_prop.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/refer_prop.py) to generate the reference IAG parameters. 
2.	Extract all the parameters / thermodynamic properties and do an averaging (if needed). 
3.	To compute isotherms, you need the natural log of partition functions (one of the outputs of the above step). You can copy them in a file. Create the “analysis” folder, copy the file to it 
4.	Run [analysis.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/analysis.py) code. One may need to modify a few things in the code:
i.	Temperature, mass of particles adsorbed, one needs to modify them as per case-by-case basis
ii.	Range of pressures/chemical potentials to investigate 
The code generates the isotherms data in the 3 files:
i.	P_predicted_from_CN.txt (prediction of pressures for considered loadings using Canonical Partition Function approach.)
ii.	N_predicted_from_GCN.txt (1st way for prediction of loading for considered pressures using Grand Canonical Partition Function approach 
iii. N_predicted_from_FL.txt ((2nd way for prediction of loading for considered pressures using Grand Canonical Partition Function approach)
Please refer "2PT Analysis of adsorption" section of our [publication](https://doi.org/10.1063/5.0099790)
5.	Plot the above to get isotherms using [comp_plots_isotherm.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/comp-plots-code.py). You may use other plotting softwares(xmgrace)
