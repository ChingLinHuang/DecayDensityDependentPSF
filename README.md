# DecayDensityDependentPSF
Data and R scripts for Huang et al.: Plant–soil feedback persists beyond host death to shape density-dependent plant competition. *Plant and Soil*. 

## Code
1. ``01_Data-Preparation.R``: Prepare ``performance.csv`` from raw data for analysis.
2. ``02_Performance_Biomass.R``: Generate Figure 2 and prepare bootstrap data for Figure 3-4 and Table S1-S6.
3. ``02_Performance_density.R``: Generate Figure S2 and prepare bootstrap data for Figure S3 and Table S7-S8.
4. ``03_plotting_Biomass.R``: Generate Figure 3-4. Prepare summary statistics for Tables S1-S6.
5. ``03_plotting_Biomass.R``: Generate Figure S3. Prepare summary statistics for Tables S7-S8.
6. ``04_tables_Biomass.R``: Generate Table S1-S6.
7. ``04_tables_density.R``: Generate Table S7-S8.
8. ``08_Performance_Biomass_different-model-fit.R``: Prepare bootstrap data for Figure S4.


## Data
**Bootstrap**
1. ``biomass-competition_boot-metrics-unnest.rdata``: Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the log(1 + competitor harvest biomass / competitor initial biomass). 
2. ``biomass-competition-08_boot-metrics-unnest.rdata``: Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the log(1 + competitor harvest biomass)
3. ``density-competition_boot-metrics-unnest.rdata``： Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the competitor density

**CleanData**
1. ``biomass_coef.rdata``: linear regression results by using log(1 + competitor harvest biomass / competitor initial biomass) as explanatory variable.
2. ``biomass-08_coef.rdata``: linear regression results by using log(1 + competitor harvest biomass) as explanatory variable.
3. ``density_coef.rdata``: linear regression results by using competitor density as explanatory variable.
4. ``performance.csv`` and ``performance.rdata``: harvest biomass measurements for all individuals in the experiment.

**Figures**
1. ``biomass-competition_comp-coef.pdf``: Figure 3b
2. ``biomass-competition_comp-outcome.pdf``: Figure 4b
3. ``biomass-competition_intrinsic-growth.pdf``: Figure 3a
4. ``biomass-competition_model-fit.pdf``: Figure 2
5. ``biomass-competition_persist.pdf``: Figure 4a
6. ``biomass-competition-08_comp-coef.pdf``: Figure S4b
7. ``biomass-competition-08_comp-outcome.pdf``: Figure S4d
8. ``biomass-competition-08_intrinsic-growth.pdf``: Figure S4a
9. ``biomass-competition-08_persist.pdf``: Figure S4c
10. ``density-competition_comp-coef.pdf``: Figure S3b
11. ``density-competition_comp-outcome.pdf``: Figure S3d
12. ``density-competition_intrinsic-growth.pdf``: Figure S3a
13. ``density-competition_model-fit.pdf``: Figure S2
14. ``density-competition_persist.pdf``: Figure S3c

**Model Selection**
1. ``biomass-competition_AIC.csv``: AIC values for models using log(1 + competitor harvest biomass / competitor initial biomass) as the explanatory variable, excluding (AIC.1) or including (AIC.2) a quadratic term.
2. ``biomass-competition_BIC.csv``: BIC values for models using log(1 + competitor harvest biomass / competitor initial biomass) as the explanatory variable, excluding (BIC.1) or including (BIC.2) a quadratic term.
3. ``biomass-competition_F-test_2_0.csv``: F-test results comparing intercept-only, linear, and quadratic models using log(1 + competitor harvest biomass / competitor initial biomass) as the explanatory variable.
4. ``density-competition_AIC.csv``: AIC values for models using competitor density as the explanatory variable, excluding (AIC.1) or including (AIC.2) a quadratic term.
5. ``density-competition_BIC.csv``: BIC values for models using competitor density as the explanatory variable, excluding (BIC.1) or including (BIC.2) a quadratic term.
6. ``density-competition_F-test_2_0.csv``: F-test results comparing intercept-only, linear, and quadratic models using competitor density as the explanatory variable.

**RawData**
1. ``harvest-dry weight.xlsx``: Harvest dry biomass of all individuals
2. ``initial biomass.xlsx``:initial dry biomass measurements for E. roxburghiana and M. zuihoensis.

**Tables**
1. ``biomass-competition_AIC_BIC.tex``: Table S2
2. ``biomass-competition_boot_a.tex``: Table S3
3. ``biomass-competition_boot_b.tex``: Table S4
4. ``biomass-competition_boot_comp.tex``: Table S6
5. ``biomass-competition_boot_persist.tex``: Table S5
6. ``biomass-competition_F-test_2_0.tex``: Table S1
7. ``density-competition_boot_comp.tex``: Table S8
8. ``density-competition_boot_persist.tex``: Table S7
