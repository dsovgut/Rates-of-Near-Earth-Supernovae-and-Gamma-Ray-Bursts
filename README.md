# Cosmic Hazard Forecast: Near-Earth Supernovae & Gamma-Ray Burst Rates 🌌

## Project Overview
This repository contains the code and datasets used in the research **"Cosmic Hazard Forecast: Rates of Near-Earth Supernovae and Gamma-Ray Bursts"** by **Brian Fields, Danylo Sovgut, Ashvini Krishnan, and Alexandra Trauth**. The study presents a model for estimating the rate **Γ(𝑟)** of cosmic explosions (supernovae & GRBs) as a function of distance from Earth.

Cosmic explosions have potential implications for Earth's biosphere, mass extinctions, and the **Galactic Habitable Zone**. The code within this repository generates datasets, performs numerical integration, and produces plots that explore these astrophysical threats.

## 🔬 Research Motivation
- **Supernovae (SNe) & GRBs** are among the most energetic cosmic events and can influence planetary environments.
- Nearby explosions may **threaten Earth's biosphere**, potentially causing **ozone depletion and increased radiation exposure**.
- Understanding the frequency and distribution of these events provides insights into **habitability conditions across the Milky Way**.

## 📊 Key Objectives
- Model the **rate of supernovae & GRBs** as a function of **distance from Earth**.
- Compare explosion rates from **different galactic distributions**:
  - **TRILEGAL model** (stellar distribution-based rate)
  - **Green's model** (SNR-based rate)
- Incorporate statistical **clustering effects** for explosion progenitors.
- Evaluate the **threat level** from cosmic explosions based on their **kill distances**.

## 🛠️ Methodology & Model Framework
### **1. Explosion Rate Density**
- The local event rate **Γ(𝑟)** is derived by integrating the **explosion rate density** over a spherical volume.
- Axisymmetric models distinguish between **thin disk (CCSN, LGRB) and thick disk (SGRB, Type Ia SN)** events.
- Rates are computed for two galaxy models:
  - **TRILEGAL Model** (follows galactic stellar density)
  - **Green Model** (based on supernova remnant distribution)

### **2. Local Explosion Rate Calculation**
- Integrates over galactocentric coordinates **(𝑅, 𝑧)** to compute the **local explosion rate near Earth**.
- Models consider **solar position and vertical oscillations in the Galactic plane**.

### **3. Clustering Effects**
- Cosmic explosions occur in **stellar clusters**, affecting the **mean time between nearby events**.
- Adjustments made to reflect clustered supernova rates.

### **4. Defining Cosmic Hazard: Kill Distances**
- A cosmic explosion is **hazardous** if it depletes **≥30% of the ozone layer**.
- Defined **kill distances** for each explosion type:
  - **Core-Collapse Supernovae (CCSN):** 10 pc
  - **Long Gamma-Ray Bursts (LGRB):** 2044 pc
  - **Short Gamma-Ray Bursts (SGRB):** 91.4 pc
  - **Type Ia Supernovae:** 10 pc

## 📂 Computational Components
### **1. Data Generation Scripts**
- `generate_trilegal_rates.py`: Computes rates for TRILEGAL model.
- `generate_isotope_rates.py`: Computes isotope model rates (Fe-60 & Al-26 contributions).
- `generate_solar_height_rates.py`: Explores CCSN & Type Ia SN rate variations with solar height.

### **2. Calculation Modules**
- `mean_times.py`: Computes mean recurrence times of explosions.
- `galactic_center_rates.py`: Estimates explosion rates near the Milky Way’s center.

### **3. Plotting Scripts**
- `rate_vs_distance_trilegal.py`: Generates explosion rate vs. distance plots.
- `rate_vs_distance_isotope.py`: Compares Trilegal vs. Isotope model.
- `rate_vs_extinctions.py`: Analyzes correlation between explosion rates and mass extinctions.
- `historical_sn_comparison.py`: Compares model predictions with historical supernova records.

## 📈 Key Findings
- The **CCSN & LGRB rates** are the most significant for potential **biospheric impact**.
- The **mean recurrence time** for a **near-Earth CCSN event** is **~2.4 Gyr**.
- A **near-Earth LGRB event** is rarer (~0.7 Gyr mean time), but far more catastrophic.
- **SGRB & Type Ia SN pose lower risks** due to their infrequent occurrence.
- **Mass extinction events may be correlated** with historical supernova activity.

## 🌍 Implications & Future Work
- The results suggest that cosmic explosions **could play a role in Earth’s extinction events**.
- This framework could be used to **map galactic habitability zones**, evaluating where planetary biospheres are least at risk.
- Future improvements:
  - Refine **stellar clustering models** for CCSN rates.
  - Integrate **spiral arm dynamics** into supernova rate distributions.
  - Extend models to include **external galaxies & cosmic ray exposure**.

---
🔹 **Authors:** *Danylo Sovgut, Brian Fields, Ashvini Krishnan, Alexandra Trauth*  
🔹 *This research is affiliated with the University of Illinois at Urbana-Champaign, Department of Astronomy.*

