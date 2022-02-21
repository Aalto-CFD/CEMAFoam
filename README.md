# CEMAFoam
![OpenFOAM 2006](https://img.shields.io/badge/OpenFOAM-2006-brightgreen) ![license](https://img.shields.io/badge/license-GPL3-brightgreen)

[CEMAFoam](https://github.com/Aalto-CFD/CEMAFoam) is a C++ library that implements Chemical Explosive Mode Analysis (CEMA) in OpenFOAM. The developments of this library are a result of the project coursework [CFD with OpenSource Software 2021](http://dx.doi.org/10.17196/OS_CFD#YEAR_2021).

## Table of contents
* [Description](#description)
* [Dependencies](#dependencies)
* [Validation](#validation)
* [How to use](#how-to-use)
* [Authors](#authors)
* [Cite](#cite)
* [References](#references)
* [License](#license)

## Description

**CEMAFoam** is a C++ library that implements Chemical Explosive Mode Analysis (CEMA) in `OpenFOAM`.
CEMA is considered a versatile computational diagnostics tool that enables the detection of various critical combustion
features including reaction fronts, flame stabilization mechanisms, and auto-ignition and extinction zones [[1]](#1).
Recently, it was further extended to account for diffusion [[2]](#2) and evaporation [[3]](#3) processes to investigate their roles toward promoting or inhibiting chemistry and auto-ignition.

The present library implements basic formulations of CEMA, which are responsible
for the identification of pre- and post-ignition zones and subsequently the reaction fronts, into
OpenFOAM. Such an implementation is considered the crucial part in the analysis tool development
since it identifies the existence of chemical explosive mode.


## Dependencies
**CEMAFoam** is currently compatible with [`OpenFOAM v2006`](https://www.openfoam.com/news/main-news/openfoam-v20-06) from ESI. Moreover, it requires [`pyJac`](https://github.com/SLACKHA/pyJac) for analytical Jacobian evaluation.

## Validation

Validation of the developed model is presented below for a laminar one-dimensional unstrained planar premixed flame comprising methane and air at equivalence ratio of 0.5 and thermodynamic conditions of T = 900 K and p = 1 atm. The adopted chemical kinetic mechanism is developed by Yao et al. [[4]](#4) comprising 54 species and 269 reactions.

The following figure depicts the flame structure through temperature and heat release fields. The variable `cem` defines the leading non-conservative eigenvalue of the thermo-chemical analytic Jacobian matrix. The positive values of `cem` indicate pre-ignition zones and the negative values indicate post-ignition zones, whereas zero-crossing interface can be regarded as the reaction front. Field plots are presented for the numerical based Jacobian (left panel) and analytical based Jacobian (right panel). 
Discrepancies are shown in the preheat zone of `cem` using numerical Jacobian are possibly due to insufficient significant digits of the Jacobian matrix resulting by finite-differencing. Such notes are further supported by discussions of the original CEMA developments by Lu et al. [[1]](#1).

![cema_valid](https://imgur.com/XgqdvPr.png)

*<div align="center">Validation of CEMA implementation for 1D methane/air laminar premixed flame. CEMA results using numerical Jacobian (left panel) show discrepancies in preheat zone indicating importance of using analytical formulation of thermo-chemical Jacobian (as indicated in right panel).</div>*

The following figure depicts local combustion modes using projections of diffusion and reaction terms onto chemical explosive mode, as proposed by Xu et al. [[2]](#2). Projected CEMA results from developed library in OpenFOAM are compared against reference implementation from PREMIX code. Implementations for the projected CEMA is not currently available in the repository but they will be uploaded soon.

![cema_valid2](https://imgur.com/yZfNicq.png)

*<div align="center">Validation of projected CEMA for combustion mode characterization. Results using OpenFOAM are compared against reference implementation using PREMIX code for the same initial conditions.</div>*


## How to use
Make sure that `OpenFOAM v2006` is installed and properly sourced, then navigate to our library directory and follow the instructions below.

* Compile the library by executing the following commands from terminal interface
```
cd src/thermophysicalModels/chemistryModel
wmake
```

* Choose the new chemistry model from the corresponding subdictionary in
  `constant/chemistryProperties` as in the following.
```
chemistryType
{
    solver            odePyjac;
    method            cemaPyjac;
}
```

* Link the library during solver runTime. This is achieved by adding the
  following to `system/controlDict` of the simulation case directory.

```
libs
(
    "libcemaPyjacChemistryModel.so"
);
```

## Authors
The open-source library is a property of [Aalto-CFD](https://github.com/Aalto-CFD) and it is developed and currently maintained by

- Mahmoud Gadalla (mahmoud.gadalla@aalto.fi)
- Islam Kabil (kabil@uconn.edu)

## Cite

If you use our model, please consider citing the following work:

<a id="cite"></a> 
M. Gadalla: Implementation of Analytical Jacobian and Chemical Explosive Mode Analysis
(CEMA) in OpenFOAM. In Proceedings of CFD with OpenSource Software, 2021, Edited by H. Nilsson. 
doi:[10.17196/OS_CFD#YEAR_2021](http://dx.doi.org/10.17196/OS_CFD#YEAR_2021).
<details>
<summary>BibTex</summary>
<p>

```
@inproceedings{Gadalla2022cema,
  author = {Gadalla, Mahmoud},
  title = {{Implementation of Analytical Jacobian and Chemical Explosive Mode Analysis (CEMA) in OpenFOAM}},
  booktitle={Proceedings of CFD with OpenSource Software},
  editor={Nilsson, Håkan},
  doi = {10.17196/OS_CFD#YEAR_2021},
  year={2021}
}
```
</p>
</details>

## References

<a id="1">[1]</a> 
T. F. Lu,  C. S. Yoo, J. H. Chen, and C. K. Law. Three-dimensional direct numerical simulation of a turbulent lifted hydrogen jet flame in heated coflow: A chemical explosive mode analysis. Journal of Fluid Mechanics, 652, 45-64. doi:[10.1017/S002211201000039X](https://doi.org/10.1017/s002211201000039x). (2010)


<a id="2">[2]</a> 
C. Xu, J.-W. Park, C. S. Yoo, J. H. Chen, and T. Lu, “Identification of premixed flame propagation modes using chemical explosive mode analysis,” Proceedings of the Combustion Institute, vol. 37, no. 2, pp. 2407–2415. doi:[10.1016/j.proci.2018.07.069](https://doi.org/10.1016/j.proci.2018.07.069). (2019)

<a id="3">[3]</a> 
D. Mohaddes, W. Xie, and M. Ihme, “Analysis of low-temperature chemistry in a turbulent swirling spray flame near lean blow-out,” Proceedings of the Combustion Institute, vol. 38, no. 2, pp. 3435–3443. doi:[10.1016/j.proci.2020.08.030](https://doi.org/10.1016/j.proci.2020.08.030). (2021)

<a id="4">[4]</a> 
T. Yao, Y. Pei, B.-J. Zhong, S. Som, T. Lu, and K. H. Luo, “A compact skeletal mechanism for n-dodecane with optimized semi-global low-temperature chemistry for diesel engine simulations,” Fuel, vol. 191, pp. 339–349. doi:[10.1016/j.fuel.2016.11.083](https://doi.org/10.1016/j.fuel.2016.11.083) (Mar. 2017)


## License

**CEMAFoam** library follows the GNU General Public License.
See the [LICENSE](LICENSE) file for license rights and limitations.
