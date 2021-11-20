# Introduction
* This work is part of the _**part decomposition (PD)**_ study [[Oh et al. (2021)](https://doi.org/10.1016/j.addma.2020.101702)].
* In the study, we developed 8 indicators to evaluate decomposed 3D models in terms of engineering design. The 8 indicators can be applied to identifying how an initial model is well decomposed into several pieces for 3D printing (also known as additive manufacturing). 
* The code helps calculate the following 8 evaluation indicators:  
  - Roughtness: the roughness of part surface
  - Overhang: the extent of critical overhang facets
  - Sharpness: the extent of shart edges
  - Gap: the area of facets causing non-allowable gaps
  - Concavity: the volume of concave spaces
  - Feasibility: the feasibility of decomposed parts
  - Interface: the area of connection interface
  - Quantity: the number of decomposed components
* The code involves two types of part decomposition (PD) algorithms: concave feature-based PD (CPD) and convex feature-based PD (VPD). When running the code an initial model is decomposed into sub-parts and then the decomposed sub-parts are evaluated based on the 8 indicators. 
  - In our study, the V-HACD algorithm is applied to CPD
  - If some of the decomposed parts are infeasible, VPD is automatically applied to make them feasible. 

![image](https://user-images.githubusercontent.com/59811007/142722693-2336bd91-45b4-499c-bd3d-2b2a865d50be.png)


![image](https://user-images.githubusercontent.com/59811007/142722692-eaa0c61d-7503-49c9-9a80-05b9822c10c3.png)

 

# Prerequisites
* [BLENDER 2.82](https://www.blender.org/download/releases/2-82/) -> **It must be Ver2.82**
* [V-HACD](https://github.com/kmammou/v-hacd)

# Usage 
(1) Install and run Blender 2.82

(2) [Import V-HACD](https://github.com/andyp123/blender_vhacd)  

(3) Put the Python code of "eval_PD.py" to the script embedded in Blender 2.82

(4) Select a 3D model 

(5) Click the "Run" button to run the script

(6) The result will be shown in Toggle System Console (MENU: Window > Toggle System Console)


# Contact
* Yosep Oh, yosep.oh@kgu.ac.kr
* Soonjo Kwon, soonjo.kwon@kumoh.ac.kr

# Version
* Version 1.0

# Reference
- Oh, Y., Ko, H., Sprock, T., Bernstein, W. Z., and Kwon, S., [**Part decomposition and evaluation based on standard design guidelines for additive manufacturability and assemblability**](https://doi.org/10.1016/j.addma.2020.101702), *Additive Manufacturing*, 37, 101702, 2021.
 
