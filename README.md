# Skin Auxetic Inverse Analysis

This repository contains MATLAB codes for inverse analysis of human back skin data and finite element simulations of skin-inspired auxetic meshing patterns.

The repository has two main purposes:

1. To fit and compare material models for human skin.
2. To simulate auxetic skin-meshing geometries using GOH and Ogden material models.

---

## Important requirement: GIBBON

These codes use the GIBBON toolbox. The simulations will not work without GIBBON because the codes use GIBBON functions for geometry generation, meshing, FEBio file creation, running FEBio analyses, and post-processing.

Please install GIBBON before running the codes:

https://www.gibboncode.org/Installation/

You also need:

- MATLAB
- FEBio
- GIBBON toolbox
- MATLAB Optimization Toolbox for the inverse analysis scripts

After installing GIBBON, make sure it is added to the MATLAB path before running the scripts.

---

## Repository structure

```text
skin-auxetic-inverse-analysis
├── inverse_analysis_of_skin_meshing
├── auxetic_geometry_GOH
└── auxetic_geometry_Ogden
```

---

## 1. Inverse analysis of skin meshing

Folder:

```text
inverse_analysis_of_skin_meshing
```

This folder contains the inverse analysis codes used to compare material models for human back skin.

The experimental stress-stretch data were extracted from Figure 10 of the human back skin study by Ní Annaidh et al. (2012):

https://link.springer.com/10.1007/s10439-012-0542-3

The data include two loading directions:

- perpendicular to the Langer line
- parallel to the Langer line

Main files:

```text
main_compare_models.m
run_GOH_2fiber.m
run_Ogden.m
```

### `main_compare_models.m`

This is the main script for model comparison.

It runs:

```matlab
res2 = run_GOH_2fiber;
res3 = run_Ogden;
```

The script then prints the model error values and plots the predicted stress-stretch curves against the extracted experimental data.

It compares:

- GOH two-fibre model
- Ogden model
- experimental perpendicular data
- experimental parallel data
- mean experimental curve

To run the inverse analysis comparison, open MATLAB, go to this folder, and run:

```matlab
main_compare_models
```

### `run_GOH_2fiber.m`

This script fits the two-fibre GOH model to the extracted human back skin data.

### `run_Ogden.m`

This script fits the Ogden model to the extracted human back skin data.

---

## 2. Auxetic geometry with GOH model

Folder:

```text
auxetic_geometry_GOH
```

This folder contains the simulation codes for auxetic skin-meshing patterns using a GOH-type material model.

Main files:

```text
run_auxetic_GOH.m
allGeometries.m
```

### `run_auxetic_GOH.m`

This is the main simulation script for the GOH auxetic geometry model.

It performs the following steps:

1. Selects one auxetic geometry.
2. Generates the geometry and mesh.
3. Applies mesh refinement.
4. Defines material directions.
5. Creates the FEBio model.
6. Runs the finite element simulation.
7. Imports and plots the results.

### `allGeometries.m`

This function stores the predefined auxetic geometry options.

The geometry is selected in the main script using:

```matlab
geometryID = 4;
geo = allGeometries(geometryID, pointSpacing, distRefine);
```

To use another geometry, change the value of `geometryID`.

For example:

```matlab
geometryID = 1;
```

or:

```matlab
geometryID = 5;
```

The available geometry options are defined inside `allGeometries.m`.

---

## 3. Auxetic geometry with Ogden model

Folder:

```text
auxetic_geometry_Ogden
```

This folder contains the corresponding auxetic geometry simulation using an Ogden material model.

Main files:

```text
run_auxetic_Ogden.m
allGeometries.m
```

The workflow is similar to the GOH folder, but the material model is different.

---

## Changing the geometry

The auxetic geometry can be changed by editing:

```matlab
geometryID = 4;
```

Possible values are:

```matlab
geometryID = 1;
geometryID = 2;
geometryID = 3;
geometryID = 4;
geometryID = 5;
```

Each value corresponds to a different predefined auxetic pattern in `allGeometries.m`.

---

## Mesh refinement and mesh convergence

Mesh density can be adjusted mainly using:

```matlab
pointSpacing = 1;
distRefine = [2];
```

### `pointSpacing`

This controls the general mesh size.

Examples:

```matlab
pointSpacing = 2;    % coarser mesh
pointSpacing = 1;    % finer mesh
pointSpacing = 0.5;  % finer mesh
```

Smaller `pointSpacing` values usually create a finer mesh, but the simulation becomes slower.

### `distRefine`

This controls local refinement near important geometric features, such as slit corners or sharp regions.

Examples:

```matlab
distRefine = [3];
distRefine = [3 2];
distRefine = [4 3 2];
```

A simple mesh convergence study can be done by running the same geometry with different `pointSpacing` and `distRefine` values, then comparing results such as force, stress, deformation, or Poisson-type response.

---

## How to run the simulations

### Inverse analysis

Go to:

```text
inverse_analysis_of_skin_meshing
```

Run:

```matlab
main_compare_models
```

### GOH auxetic simulation

Go to:

```text
auxetic_geometry_GOH
```

Run:

```matlab
run_auxetic_GOH
```

### Ogden auxetic simulation

Go to:

```text
auxetic_geometry_Ogden
```

Run:

```matlab
run_auxetic_Ogden
```

---

## Output files

The scripts may generate temporary FEBio input files, log files, and simulation output files.

These files can usually be regenerated by rerunning the MATLAB scripts, so large temporary output files are not included in this repository.

---

## Notes for users

Before running the codes, check that:

1. GIBBON is installed.
2. FEBio is installed.
3. GIBBON is added to the MATLAB path.
4. FEBio can be called from MATLAB.
5. The required `.m` files are in the same folder or added to the MATLAB path.

---

## Research context

This repository supports computational studies of skin meshing and auxetic pattern behaviour. The codes are intended to help study how geometry and constitutive modelling affect deformation, stress response, and expansion behaviour in skin-inspired structures.

---

## Author

Masoumeh Razaghi  
PhD candidate in computational biomechanics  
Atlantic Technological University, Galway
