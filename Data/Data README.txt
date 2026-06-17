--- Bootstrap
biomass-competition_boot-metrics-unnest.rdata: Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the log(1 + competitor harvest biomass / competitor initial biomass)

biomass-competition-08_boot-metrics-unnest.rdata: Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the log(1 + competitor harvest biomass)

density-competition_boot-metrics-unnest.rdata： Bootstrap samples of single individual growth, and competitor effects when regressing focal species growth on the competitor density　

--- CleanData
biomass_coef.rdata: linear regression results by using log(1 + competitor harvest biomass / competitor initial biomass) as explanatory variable.

biomass-08_coef.rdata: linear regression results by using log(1 + competitor harvest biomass) as explanatory variable.

density_coef.rdata: linear regression results by using competitor density as explanatory variable.

performance.csv and performance.rdata: harvest biomass measurements for all individuals in the experiment.

--- Figures

Contains all figures used in the manuscript.

--- Model Selection
biomass- or density-: using log(1 + competitor harvest biomass / competitor initial biomass) or competitor density as explanatory variable.

_AIC: AIC values when the model excludes (AIC.1) or include (AIC.2) quadratic form of explanatory variable.

_BIC: BIC values when the model excludes (BIC.1) or include (BIC.2) quadratic form of explanatory variable.

_F-test: F-test results comparing intercept-only, linear, and quadratic models.

--- RawData
Harvest dry biomass of all individuals and initial dry biomass measurements for E. roxburghiana and M. zuihoensis.

--- Tables
biomass- or density-: using log(1 + competitor harvest biomass / competitor initial biomass) or competitor density as explanatory variable.

-AIC_BIC and F-test2_0: Results of model selection 

-boot_a: bootstrap distributions of single individual growth

-boot_b: bootstrap distributions of competitor effects

-boot_comp: bootstrap distributions of competitive outcomes

-boot_persist: bootstrap distributions of persistence


--- 