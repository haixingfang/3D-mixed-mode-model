# 3D-mixed-mode-model
The codes were developed for implementing a 3D mixed-mode model to predict the average ferrite grain size and grain size distribution for an austenite-to-ferrite phase transformation during continuous cooling of an Fe-C-Mn steel. Using a Voronoi construction to represent the austenite grains, the ferrite is assumed to nucleate at the grain corners and to grow as spheres. Classical nucleation theory is used to estimate the density of ferrite nuclei. By assuming a negligible partition of manganese, the moving ferrite–austenite interface is treated with a mixed-mode model in which the soft impingement of the carbon diffusion fields is considered. The ferrite volume fraction, the average ferrite grain size, and the ferrite grain size distribution are derived as a function of temperature. The model provides a versatile tool to analyze the evolution of the ferrite grain size distribution at low computational costs.

# Features of the model
- Continuous nucleation at the corners of austenitic grains
- Interface moving in a mixed-mode manner by assuming a negligible partition of manganese in para-equilibrium
- Using effective interface mobility to account for the influence of substitutional alloying elements
- Soft impingement and hard impingement are considered

# Denpendencies of the code
- Install the Multi-Parametric Toolbox 3 (mpt3): https://www.mpt3.org/ into the same folder, for generating Voronoi cells to represent austenite grains.
- Make sure that 'Optimization toolbox' and 'Symbolic Math toolbox' are included in your own Matlab package. These toolboxes should be included by default.
All codes have been tested executable with Matlab 2014b or above.

# How to run the code
Just run [ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m). <br>
But, remember to first set up the chemical compositions, heat treatment parameters and include thermodynamic data (which can be parameterized first with [Thermo-Calc software](https://www.thermocalc.com/)) in [SimulCond.m](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/SimulCond.m). <br>
<br>

# License
This package is free to use, ditribute and adapt, but no warranty and liability to any kinds of simulation results. <br>
Citing our [article](https://link.springer.com/content/pdf/10.1007/s11661-017-4397-y.pdf) is strongly encouraged if you use or get inspired by our code. <br>
See the [LICENSE](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/LICENSE) for license rights and limitations (GNU General Public License v3.0).

## Contact via hfang@tudelft.nl or haixingfang868@gmail.com
