# Rates-of-Near-Earth-Supernovae-and-Gamma-Ray-Bursts
This code repository includes all the coding used to generate results in the paper Cosmic Hazard Forecast:  Rates of Near-Earth Supernovae and Gamma-Ray Bursts. 
The codes included generate datasets, plots, and calculations displayed in the paper. 
At the top of the each code, the inputs such as rates for cosmic explosions or scale heights and radius are defined. 
If one wishes to use our model with different inputs, he/she only needs to change the inputs and run the code to generate
the plots and results with his/her inputs. 
### Prerequisites
All the codes are generated in Python and no other software is needed in addition to standart 
Python packages (such as matplotlib or numpy that already come with Anaconda) are needed.

## 1. Generating Datasets for Plots:
This folder includes three files that generate rates for cosmic explosions at various distances. 
Since the full integration of the rate at a given distance is computationally intensive and time consuming, 
it makes sense to make the plots from the datasets directly. Each of these files generates a dataset for our plots. 
However, if one wishes to change the inputs to our model and generate and new dataset and make a new plot, he/she can use 
these codes. <br />
Each file generates two csv datasets. One file contains rates for the thick and thin height and radius. In the appendix of our paper,
we describe in-depth why that is the only calculation needed to generate the cosmic explosion rates.
Each cosmic explosion is comprised of the thick or thin component of the rate multiplied by the rate of the given explosion.
Generating the dataset this way significantly speeds up the calculation since we only need to generate the two components 
instead of the 4 different rates for cosmic explosions. Core-Collapse Supernovae (CCSN)
and Long Gamma-Ray Bursts (GRB) follow the thin disk, Short Gamma-Ray Bursts (SGRB) follows thick disk, and 
Type Ia supernovae are equally split between the thick and thin disk component.
Here is an explanation of what type of dataset each file generates:
### a) Generating Dataset [Trilegal Rates]
This file generates scale thick and thin disk components for the Trilegal Model. It also converts it to the rate for 
each type of explosion at a given distance.  
### b)Generating Dataset [Solar Heights Rates]
This file generates a dataset from the Trilegal Model of the CCSN and Type Ia. The rate for each explosion is generated
varying with the solar height. 
### c) Generating Dataset [Isotope Rates] 
This file generates scale thick and thin disk components for the Isotope Model both Al26 and Fe60. It also converts it to the rate for 
each type of explosion at a given distance.  
## Datasets
In order for all the codes to work, each csv file in the directory needs to be in the same folder
as the code or os.chdir function needs to be added to the codes making plots if the directory of csv files
is different from the directory of csv files on your local machine. Below is a quick explanation of what 
each csv files contains:
1. scale_rates2.csv - scale rates for the Trilegal model
<br />
2. final_rates2.csv - rates for each type of cosmic explosion for the Trilegal Model
<br />
3. Scale_rates_aluminum_ferum.csv - scale rates for the isotope model
<br />
4. rates-Al-Fe.csv - rates for each type of cosmic explosion for the Isotope Model
<br />
5. CCSN_zheights.csv - rates for CCSN for different solar heights for the Trilegal Model
<br />
6. SN_zheights -> Type Ia rates for different solar heights for the Trilegal model

## Conducting Calculations
These files conduct calculations mentioned in our tables and results section. No csv files needed
to replicate the results. 
### a)Calculating Mean Times for Cosmic Explosions - calculates mean times with errors for each type of 
cosmic explosion. Results from this file are referenced in Table 4 of our paper. Calculations are made
for the Trilegal Model, Isotope Model, and the SNR Model
### b)Calculating Rate at the Galactic Center - calculates the rate of each type of explosion at the galactic center.
This result is referenced in our results section. Particularly useful if one wishes to quickly calculate the rate 
at a given location in the Milky Way. 
## Making Plots 
Each of these codes makes a plot that is referenced in our paper. Most of the use the afformentioned datasets. 
### a) Rate vs Distance Plot [TRILEGAL, Approximations]
Makes a plot of rate vs distance for the cosmic explosions for the Trilegal Model. Also compares
our approximation functions to the full integration. If anyone is interested in the approximation functions which 
we discuss in detail in the appendix of our paper, this is the file for these functions.
### b) Rate vs Distance Plot [Isotope vs Trilegal Model]
The rate vs distance plot comparing the Trilegal Model with the Isotope Model.
### c) Rate vs Distance Plot [Solar Heights]
Rate vs Distance plot for CCSN and Type Ia supernovae with different solar heights. 
### d) Rate vs Mass Exctinctions Plot
Plots the rate of cosmic explosions vs the rate of mass exctinctions. 
### e) Rate vs Distance Plot [Historical Supernovae]
Compares the rate vs distance for the TRILEGAL model with the rate of known observed Historical Supernovae.

## Authors
Danylo Sovgut and Brian Fields
