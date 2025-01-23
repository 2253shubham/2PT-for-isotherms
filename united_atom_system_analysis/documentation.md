# Quick Steps to Perform 2PT Analysis

## 1. Generate Velocity Autocorrelation Data
- Perform MD simulations to generate velocity autocorrelation data.
- Refer to the **"Computational details for methane adsorption in silicalite"** section of our [publication](https://doi.org/10.1063/5.0099790) to understand how to set up the simulations and required parameters.

---

## 2. Run Density of States (DoS) Computation Scripts
- Use the velocity autocorrelation data obtained from MD simulations and run one of the following scripts:
  - [IAG_DoS_comp.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/HS_DoS_comp.py) (for Ideal Adsorbed Gas systems)
  - [HS_DoS_comp.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Hard_Sphere_approximation/scripts/HS_DoS_comp.py) (for Hard Sphere systems)

### Note for IAG Systems:
- Run the additional script [refer_prop.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/refer_prop.py) to generate reference IAG parameters.

---

## 3. Extract Parameters and Thermodynamic Properties
- Extract all parameters and thermodynamic properties.
- Perform averaging if required.

---

## 4. Compute Isotherms
1. Obtain the natural log of partition functions (one of the outputs from the previous step).
2. Save the results in a file and create an `analysis` folder. Copy the file to this folder.
3. Run the [analysis.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/analysis.py) script. 

### Modifications Needed:
- Update the script with case-specific parameters:
  - **Temperature**
  - **Mass of adsorbed particles**
  - **Range of pressures/chemical potentials**

### Output Files:
- The script generates the following files with isotherms data:
  1. `P_predicted_from_CN.txt`: Prediction of pressures for considered loadings using the Canonical Partition Function (CN) approach.
  2. `N_predicted_from_GCN.txt`: Loading predictions for considered pressures using the Grand Canonical Partition Function (GCN) approach (Method 1).
  3. `N_predicted_from_FL.txt`: Loading predictions for considered pressures using the GCN Function approach (Method 2).

Refer to the **"2PT Analysis of adsorption"** section of our [publication](https://doi.org/10.1063/5.0099790) for further details.

---

## 5. Plot Isotherms
- Use the output files to generate isotherm plots with the script [comp_plots_isotherm.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/united_atom_system_analysis/Ideal_Adsorbed_Gas_approximation/scripts/comp-plots-code.py).
- Alternatively, other plotting software like **xmgrace** can also be used.

---
