# [3D-mixed-mode-model](https://github.com/haixingfang/3D-mixed-mode-model)
## Austenite-ferrite TransModel 1.0
The codes were developed for implementing a 3D mixed-mode model to predict the average ferrite grain size and grain size distribution for an austenite-to-ferrite phase transformation during continuous cooling of an Fe-C-Mn steel. Using a Voronoi construction to represent the austenite grains, the ferrite is assumed to nucleate at the grain corners and to grow as spheres. Classical nucleation theory is used to estimate the density of ferrite nuclei. By assuming a negligible partition of manganese, the moving ferrite–austenite interface is treated with a mixed-mode model in which the soft impingement of the carbon diffusion fields is considered. The ferrite volume fraction, the average ferrite grain size, and the ferrite grain size distribution are derived as a function of temperature. The model provides a versatile tool to analyze the evolution of the ferrite grain size distribution at low computational costs. The codes were a result of [Haixing Fang](https://orcid.org/0000-0001-8114-5276)'s [PhD thesis](https://repository.tudelft.nl/islandora/object/uuid%3Aecd8e101-3164-4227-b47b-13a04bc4b8fb?collection=research) supervised by [Dr.ir. N.H. van Dijk](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/radiation-science-technology/research/research-groups/fundamental-aspects-of-materials-and-energy/people/niels-van-dijk/) and [Prof.dr.ir. S. van der Zwaag](https://www.tudelft.nl/lr/organisatie/afdelingen/aerospace-structures-and-materials/novel-aerospace-materials/people/personal-pages-novam/academic-staff/s-van-der-zwaag-sybrand/) at Delft University of Technology.

# Features of the model
- Nucleation at the corners of austenitic grains based on continuous nucleation theory (CNT) or simplified nucleation model (SNM) 
- Interface moving in a mixed-mode manner by assuming a negligible partition of manganese in para-equilibrium
- The influence of substitutional alloying elements is accounted by an approching using effective interface mobility 
- Soft impingement and hard impingement are considered

# Denpendencies of the code
- Install the Multi-Parametric Toolbox 3 (mpt3): https://www.mpt3.org/ into the same folder, for generating Voronoi cells to represent austenite grains.
- Make sure that 'Optimization toolbox' and 'Symbolic Math toolbox' are included in your own Matlab package. These toolboxes should be included by default.
All codes have been tested executable with Matlab 2014b or above.

# How to run the code
Just run [ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m](https://github.com/haixingfang/3D-mixed-mode-model/blob/master/ferrite_3d_model_voronoin_PBC_ND_CNT.m). <br>
Currently, it is only valid for simulating phase transformations during continuous cooling for Fe-0.1C-0.49Mn (wt%). <br>
But it can be easily adapted to other alloys under different heat treatment conditions.

# License
This package is free to use, ditribute and adapt, but no warranty and liability to any kinds of simulation results. <br>
See the [LICENSE](https://github.com/haixingfang/3D-mixed-mode-model/blob/master/LICENSE) for license rights and limitations (GNU General Public License v3.0).

# Reference
H. Fang, M.G. Mecozzi, E. Brϋck, S. van der Zwaag, N.H. van Dijk, Analysis of the grain size evolution for ferrite formation in Fe-C-Mn steels using a 3D model under a mixed-mode interface condition, Metall. Mater. Trans. A 49 (2018) 41-53. <br>
Citing our [article published in MMTA](https://link.springer.com/content/pdf/10.1007/s11661-017-4397-y.pdf) is strongly encouraged if you use or get inspired by our code.

## Contact via hfang@tudelft.nl or haixingfang868@gmail.com
