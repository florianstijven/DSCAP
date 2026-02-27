# Number of bootstrap replications and replications for the permutation test.
B = 20
# Number of MC replications for the simulations 
N_MC = 10
# Number of MC IPD replicates (per trial) to approximate true values.
n_MC <- 1e5
# Rscripts that contain helper functions.
analysishelpers = R/helper-functions/DSCAP-estimators.R R/helper-functions/permutation-LRT.R R/helper-functions/treatment-effect-estimators.R
simulationhelpers = R/simulations/generate_simulated_data.R R/simulations/parameter_values.R 
# Filename of the processed data. This can either be the original processed data
# or the synthetic processed data. 
data = data/processed_data.csv
# Formula for the logistic and linear regression models used for standardizing 
#treatment effects to the target trial.
formula = risk_score+age.geq.65+riskxage+BMI_underweight_normal

.PHONY: all application simulation

all: application simulation

simulation: results/simulations/raw-results/simulations_results_tbl.rds \
	results/simulations/raw-results/true_values_tbl.rds

application: R/neut_AZ_full_M1_estwts.Rout R/spike_AZ_full_M1_estwts.Rout \
	R/neut_AZ_truncated_M2_estwts.Rout R/spike_AZ_truncated_M2_estwts.Rout \
	R/plots_tables.Rout
# R/data-exploration.Rout
# R/data-exploration.R can only be run when the original data are available. We 
# therefore commented out this target. If the original data are available, then 
# one can uncomment this target and run make. 

# R/data-exploration.Rout: $(data)
# 	Rscript R/data-exploration.R > $@ 2> $@


# Analyses with all trials. 

R/neut_AZ_full_M1_estwts.Rout: R/estimate_dscap.R $(analysishelpers) $(data)
	Rscript --verbose R/estimate_dscap.R neut AZ $(formula) 1 $(B) $(B) 1 0 $(data) > $@ 2> $@
	
R/spike_AZ_full_M1_estwts.Rout: R/estimate_dscap.R $(analysishelpers) $(data)
	Rscript --verbose R/estimate_dscap.R spike AZ $(formula) 1 $(B) $(B) 1 0 $(data) > $@ 2> $@
	
# Analyses with J&J (Colombia), J&J (Brazil), and Novavax left out 

R/neut_AZ_truncated_M2_estwts.Rout: R/estimate_dscap.R $(analysishelpers) $(data)
	Rscript --verbose R/estimate_dscap.R neut AZ $(formula) 2 $(B) $(B) 1 1 $(data) > $@ 2> $@
	
R/spike_AZ_truncated_M2_estwts.Rout: R/estimate_dscap.R $(analysishelpers) $(data)
	Rscript --verbose R/estimate_dscap.R spike AZ $(formula) 2 $(B) $(B) 1 1 $(data) > $@ 2> $@


# Generate all plots and tables that summarize the results of the analyses. 
R/plots_tables.Rout: R/neut_AZ_full_M1_estwts.Rout R/spike_AZ_full_M1_estwts.Rout \
	R/neut_AZ_truncated_M2_estwts.Rout R/spike_AZ_truncated_M2_estwts.Rout
	Rscript --verbose R/plots_tables.R AZ > $@ 2> $@
	
	
results/simulations/raw-results/simulations_results_tbl.rds: R/simulations/run_simulation.R $(analysishelpers) $(simulationhelpers)
	Rscript R/simulations/run_simulation.R $(N_MC) $(B) $(B)
	
results/simulations/raw-results/true_values_tbl.rds: $(analysishelpers) $(simulationhelpers)
	Rscript results/simulations/raw-results/true_values_tbl.rds $(n_MC)
