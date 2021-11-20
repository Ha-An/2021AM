# Introduction
* This work is part of the _**part decomposity (PD)**_ study [[Oh et al. (2021)](https://doi.org/10.1016/j.addma.2020.101702)].
* In the study, we developed 8 indicators to evaluate decomposed 3D models in terms of engineering design.
* The 8 indicators can be applied to identifying how an initial model is well decomposed into several pieces for 3D printing (also known as additive manufacturing). 
* The code helps calculate the following 8 evaluation indicators:  
  - Roughtness: the roughness of part surface
  - Overhang: the extent of critical overhang facets
  - Sharpness: the extent of shart edges
  - Gap: the area of facets causing non-allowable gaps
  - Concavity: the volume of concave spaces
  - Feasibility: the feasibility of decomposed parts
  - Interface: the area of connection interface
  - Quantity: the number of decomposed components


# Prerequisites
* [BLENDER 2.82](https://www.blender.org/download/releases/2-82/) -> **It must be Ver2.82**
* [V-HACD](https://github.com/kmammou/v-hacd)

# Usage
* The Python code of "Evaluation_indicators" is run in the script embedded in Blender 2.82

# Contact Information
* Yosep Oh, yosep.oh@kyonggi.ac.kr
* Soonjo Kwon, soonjo.kwon@kumoh.ac.kr

# Version
* Version 1.0

# Reference
- Oh, Y., Ko, H., Sprock, T., Bernstein, W. Z., and Kwon, S., [**Part decomposition and evaluation based on standard design guidelines for additive manufacturability and assemblability**](https://doi.org/10.1016/j.addma.2020.101702), *Additive Manufacturing*, 37, 101702, 2021.
 
