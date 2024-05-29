# Amphibole compositions record cold post-emplacement reequilibration in plutons

## Data

Compiled and simulated data used in this study are included.

### Description of the data and file structure

- ds01_Whole_Rock_and_Temperatures_Pressures_Equilibrium_Melt_Compositions_Higgins2021.csv.zip

Includes whole-rock compositional data for magmatic rocks (volcanic and intrusive) compiled from the GEOROC database and other published data, chemical compositional data for hornblende, and temperature, pressure, and equilibrium melt compositions calculated by the method in Higgins et al. (2021).


- ds02_Whole_Rock_and_Temperatures_Putirka2016_Pressures_Ridolfi2022_Equilibrium_Melt_Compositions_Zhang2017.csv.zip

Includes whole-rock compositional data for magmatic rocks (volcanic and intrusive) compiled from the GEOROC database and other published data, chemical compositional data for hornblende, and temperature (Putirka2016), pressure (Ridolfi2022), and equilibrium melt compositions calculated by the method in Zhang et al. (2017).


- ds03_Whole_Rock_Composition_and_Temperatures_Pressures_Derived_from_Biotite_LiandZhang.csv.zip

Includes whole-rock compositional data for magmatic rocks (volcanic and intrusive) compiled from the GEOROC database and other published data, chemical compositional data for biotite, and temperature and pressure calculated by the method in Li & Zhang (2022).


- ds04_Simulation_Results_Using_Big_Dataset_as_Initial_Composition.csv.zip

Simulated dataset using whole-rock compositions of GEOROC arc magmatic rocks as initial compositions calculated by Perple_X, including melt, solid, and whole-rock chemical compositions.

### Sharing/Access information

Data was derived from the following sources:

- GEOROC: [https://georoc.mpch-mainz.gwdg.de//georoc/](https://georoc.mpch-mainz.gwdg.de//georoc/)

## Code

- bulk_hbl_melt_Higgins.py

Comparison of whole-rock compositions and melt compositions in equilibrium with hornblende in natural volcanic and plutonic rocks.

- amph_bulk_mpi_georoc_vol_all.jl

Used to simulate the magma evolution process on a high-performance computer cluster, generating the simulation data for this study.

- hbl_melt_bulk_compare.py

Comparison of the chemical behavior of simulated volcanic and plutonic rocks in cooling crystallization, melt extraction, and reequilibration modes.

- silica-temperature-modeling.py

Comparison of natural volcanic and plutonic rock temperatures calculated from hornblende thermometers and simulated volcanic and plutonic rock temperatures in cooling crystallization, melt extraction, and reequilibration modes as a function of silica.
