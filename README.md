# 2PT to compute adsorption isotherms

This repository contains the scripts to perform 2PT (Two Phase Thermodynamic) calculations to compute thermodynamic properties and isotherms from molecular dynamics (MD) simulations of gas-zeolite adsorption systems. 

The [2PT method](https://doi.org/10.1063/1.1624057) was developed for bulk fluid phases treats the gas-like components as hard spheres (HSs), which correctly recovers the limiting behaviors of unconfined fluids. We showed that this treatment, however, does not always lead to the correct zero-loading behavior in strongly confining systems. For methane adsorption into zeolite MFI, the HS reference state underestimates entropy by up to 20% at low loadings and leads to an order-of-magnitude increase in the adsorption onset pressure. To fix these issues, we propose the use of ideal adsorbed gas (IAG) as the gas reference model, the properties of which can be computed using the Widom insertion method on an empty adsorbent. You can read all about it in our [publication](https://doi.org/10.1063/5.0099790).

Here, we provide scripts to perform 2PT analysis on 2 systems: adsorption of molecules where molecules are modelled as single interaction sites using United Atom Model (UAM) ([united_atom_system_analysis](https://github.com/2253shubham/2PT-for-isotherms/tree/dev/united_atom_system_analysis)), and adsorption of molecules, modelled as they are ([molecular_system_analysis](https://github.com/2253shubham/2PT-for-isotherms/tree/dev/molecular_system_analysis)). 


## Documentation

[Documentation](https://github.com/2253shubham/2PT-for-isotherms/tree/dev/united_atom_system_analysis/documentation.md) for [united_atom_system_analysis](https://github.com/2253shubham/2PT-for-isotherms/tree/dev/united_atom_system_analysis)


[Documentation](https://github.com/2253shubham/2PT-for-isotherms/blob/dev/molecular_system_analysis/documentation.md) for [molecular_system_analysis](https://github.com/2253shubham/2PT-for-isotherms/tree/dev/molecular_system_analysis)


## Authors

- [@shubham](https://github.com/2253shubham)


## Environment Variables

Required python version  - 3.9 or higher \
Download and install [requirements.txt](https://github.com/2253shubham/2PT-for-isotherms/blob/main/requirements.txt) to install python library dependencies.
```bash
pip install -r /path/to/requirements.txt
```

## Appendix

If you are using this work, please cite: 

[S. Malviya, J.C. Tapia, and P. Bai, “Calculating adsorption isotherms using the two-phase thermodynamic method and molecular dynamics simulations,” J. Appl. Phys. 132(3), 034701 (2022)](https://doi.org/10.1063/5.0099790).
