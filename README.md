# Particle Image Cross-Correlation Spectroscopy (PICCS)
==========

PICCS is a spatial analysis algorithm implemented in **MATLAB**, designed specifically for Single-Molecule Localization Microscopy (SMLM) data. It quantifies the degree of co-localization or clustering between two molecular species by calculating and fitting the **Cumulative Cross-Correlation Function** according to the publication by Semrau et al. (2011).

---

## Package Components

This package consists of the following interdependent MATLAB functions, forming a complete analysis and validation pipeline:

| File Name | Purpose |
| :--- | :--- |
| `simulate_complex_formation...m` | **Simulation Engine.** Generates synthetic $A$ and $B$ localizations with user-defined densities, localization precision and percentage of binary complexes (interaction strength). |
| `PICCS_calculate_cum_corr.m` | **Calculation Engine.** Computes the cumulative cross-correlation ($C_{cum}$) between two sets of spatial coordinates. |
| `PICCS_uniform_bckgrnd.m` | **Fitting Engine.** Fits the theoretical model to $C_{cum}$ to estimate interaction strength, background particle density and localization precision. |

---

## Installation and Setup

### Prerequisites

1.  **MATLAB** installation (tested with 2024a).
2.  MATLAB's **Statistics and Machine Learning Toolbox**.
3.  MATLAB's **Optimization Toolbox**.

### Cloning the Repository

To get started, clone the repository to your local machine:

```bash
git clone git@github.com:christianpaolorichter/PICCS.git
```

or using direct download:

1.  **Download ZIP:** Navigate to the repository page on GitHub. Click the **green `<> Code`** button and select **"Download ZIP."**
2.  **Unzip:** Extract the contents of the downloaded ZIP file to your desired project location (e.g., your MATLAB projects folder).

Then, add the main directory and all subfolders to your MATLAB environment:

```Matlab
% Run this command in MATLAB after cloning:
addpath(genpath('/path/to/PICCS/'));
```

## Getting Started

To execute the toydata simulation and subsequent PICCS analysis:

```Matlab
% Reset the global random number stream to ensure reproducibility of the simulation.
clear all
reset(RandStream.getGlobalStream)

% --- Simulation Parameters ---
numCell = 1000; % total number observed cells (Number of Monte Carlo trials/simulations)
area = 	20^2;   %[µm^2] typical observable cell area
molDens = [0.4 0.35]; %[µm^-2] A >= B (Molecular density for species A and B)
prcntComplex = 0.035; %[x100% B] (Target fraction of B molecules forming complexes)
% NOTE: The final, actual percentage of complexes may slightly deviate from this target
% because the calculated number of complexes (nAB) must be rounded to a discrete integer.
epsilon = [0.025 0.025]; %[µm] (Localization precision (epsilon) for species A and B)
corrLimit = 0.2; %[µm] (Maximum radius/distance for correlation analysis)
dr = 0.01; %[µm] (Radial step size/bin width for calculating Ccum)
verbose = 1; Boolean flag (1 or 0) to display the fitted Ccum

% --- Simulation Loop ---
% Loop runs 'numCell' times, generating independent localization data sets.
for idxCell = numCell:-1:1
    % Generate synthetic localization data (A and B coordinates).
    % The function simulates a uniform distribution with a defined percentage of binary complexes.
    [A,B,prcntComplex] = ...
        simulate_complex_formation_uniform_background(...
        area,molDens,prcntComplex,epsilon,corrLimit);

    %% [AB]/[B] & [A]

    % Define the vector of radial distances (r) at which the CCF will be calculated.
    rCorr = dr:dr:corrLimit;

    % Calculate the Cumulative Cross-Correlation for the current cell.
    Ccum(:,idxCell) = PICCS_calculate_cum_corr(A,B,rCorr);
end %fun

% --- Analysis of Averaged Data ---
% The parameters are estimated by fitting the theoretical model to the mean CCF profile.
PICCS_uniform_bckgrnd(mean(Ccum,2),rCorr,verbose);
```

You should obtain the following Output demonstrating successful parameter recovery:

<p align="center">
  <img src="https://github.com/christianpaolorichter/PICCS/blob/main/results_evaluate_PICCS.png?raw=true" alt=""/>
</p>

## Citation
If you use this software/repository for your research, please cite:

### 1. The Original Publication
Semrau S, Holtzer L, González-Gaitán M, Schmidt T. 
Quantification of biological interactions with particle image cross-correlation spectroscopy (PICCS). 
Biophys J. 2011 Apr 6;100(7):1810-8. doi: 10.1016/j.bpj.2010.12.3746. PMID: 21463595; PMCID: PMC3072609.

### 2. The PICCS Software (v1.0.0)
Richter, C.P (2025). *PICCS* (v1.0.0). [Software]. Available from: https://github.com/christianpaolorichter/PICCS