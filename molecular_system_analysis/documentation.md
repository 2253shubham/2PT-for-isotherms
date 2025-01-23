# Quick Steps to Perform 2PT Analysis

> **Note**: All thermodynamic properties hardcoded in the scripts pertain to water treated using the TIP4P model, as it was the adsorbate used.

## 1. Methods to Compute Density of States (DoS)

### 1.1 Using FFT/IFFT (Indirect Method)
- **Code Location**: [Scripts for FFT/IFFT](https://github.com/2253shubham/2PT-for-isotherms/tree/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts)
  
**Steps:**
1. Compute the Velocity Autocorrelation Function (VACF) for each atom:
   - Perform FFT of the velocity array (pad with zeros first).
   - Take the modulus, square the result, and compute the inverse FFT.
   - Consider the real part of the inverse FFT as the VACF (length = `2n`, where `n` is the length of the velocity array).
2. Extract the first `n/2` correlation, mirror it, and compute the DoS using FFT on the mirrored VACF.
3. Multiply the VACF of each atom by its mass and compute the total VACF (sum of all atomic VACFs). DoS is computed per molecule.

---

### 1.2 Without FFT/IFFT (Direct Method)
- **Code Location**: [Scripts for Direct Method](https://github.com/2253shubham/2PT-for-isotherms/tree/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts)
  
**Steps:**
1. Perform FFT of the velocity array for each atom.
2. Compute the modulus, square the result, and normalize by the array length.
3. Multiply the DoS by the corresponding atomic masses.

---

## 2. Performing 2PT Analysis

### 2.1 Indirect Method (FFT/IFFT)
#### Level of Treatment for DoS Components:
1. **Translational**: 
   - Solid: Classical Harmonic Oscillator (CHO)
   - Gas: Ideal Adsorbed Gas (IAG)
2. **Rotational**: 
   - Solid: CHO
   - Gas: Hard Sphere (HS)

#### Steps:
1. Run the main code: [dos_H2O.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/dos_H2O.py).
   - A sample command is provided at the top of the file.
   - To loop over multiple systems, execute the script: [run_script.sh](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/run_script.sh).
2. Compute reference state properties using: [dos_H2O_ref.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/dos_H2O_ref.py).
   - Use the script: [ref_script.sh](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/ref_script.sh).
3. Extract thermodynamic properties (e.g., natural log of partition functions) from the outputs.
4. Compute isotherms using: [analysis.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/analysis.py).

---

### 2.2 Direct Method (Without FFT/IFFT)
#### Level of Treatment for DoS Components:
1. **Translational**:
   - Solid: Classical/Quantum Harmonic Oscillator (CHO/QHO)
   - Gas: IAG
2. **Rotational**:
   - Solid: CHO/QHO
   - Gas: HS

#### Steps:
1. Run the main code: [dos_computation.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/dos_computation.py).
   - Execute the script: [dos.computation.sh](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/dos_computation.sh).
2. **For Translational Gas Reference (IAG)**:
   - Compute reference properties using: [ref_therm_prop_calc.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/Reference_Ideal_Adsorbed_Gas/ref_therm_prop_calc.py).
   - Compute system properties using: [therm2_prop_calc.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/Reference_Ideal_Adsorbed_Gas/therm2_prop_calc.py).
   - Execute both using: [run_script_2.sh](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/Reference_Ideal_Adsorbed_Gas/run_script_2.sh).
3. **For Translational Gas Reference (HS)**:
   - Compute system properties using: [therm2_prop_calc.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/Reference_Hard_Sphere/therm2_prop_calc.py).
   - Execute using: [run_script.sh](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/without_FFT_IFFT/scripts/Reference_Hard_Sphere/run_script.sh).
4. Extract thermodynamic properties (e.g., natural log of partition functions assuming CHO and QHO).
5. Compute isotherms using: [analysis.py](https://github.com/2253shubham/2PT-for-isotherms/blob/main/molecular_system_analysis/water_in_zeolites/with_FFT_IFFT/scripts/analysis.py).

---

## 3. Modifications Required
- Input parameters in `analysis.py`:
  - Temperature
  - Mass of particles adsorbed
  - Number of molecules adsorbed
  - File locations for data
- Adjust these parameters as per the specific system being analyzed.
