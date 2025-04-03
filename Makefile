# Number of bootstrap replications and replications for the permutation test.
B = 2000
# Rscripts that contain helper functions.
helpers = R/helper-functions/estFUN.R R/helper-functions/runDSCAP.R
# Filename of the processed data. This can either be the original processed data
# or the synthetic processed data. 
data = data/processed_data_synthetic.csv
# Formula for the logistic and linear regression models used for standardizing 
#treatment effects to the target trial.
formula = risk_score+age.geq.65+riskxage+BMI_underweight_normal

.PHONY: all
all: R/neut_AZ_full_M1_estwts.Rout R/spike_AZ_full_M1_estwts.Rout \
	R/neut_AZ_truncated_M2_estwts.Rout R/spike_AZ_truncated_M2_estwts.Rout \
	R/plots_tables.Rout
# R/data-exploration.Rout
	
# R/data-exploration.R can only be run when the original data are available. We 
# therefore commented out this target. If the original data are available, then 
# one can uncomment this target and run make. 

# R/data-exploration.Rout: $(data)
# 	Rscript R/data-exploration.R > $@ 2> $@


# Analyses with all trials. 

R/neut_AZ_full_M1_estwts.Rout: $(helpers) $(data)
	Rscript --verbose R/estimate_dscap.R neut AZ $(formula) 1 $(B) $(B) 1 0 $(data) > $@ 2> $@
	
R/spike_AZ_full_M1_estwts.Rout: $(helpers) $(data)
	Rscript --verbose R/estimate_dscap.R spike AZ $(formula) 1 $(B) $(B) 1 0 $(data) > $@ 2> $@
	
# Analyses with J&J (Colombia), J&J (Brazil), and Novavax left out 

R/neut_AZ_truncated_M2_estwts.Rout: $(helpers) $(data)
	Rscript --verbose R/estimate_dscap.R neut AZ $(formula) 2 $(B) $(B) 1 1 $(data) > $@ 2> $@
	
R/spike_AZ_truncated_M2_estwts.Rout: $(helpers) $(data)
	Rscript --verbose R/estimate_dscap.R spike AZ $(formula) 2 $(B) $(B) 1 1 $(data) > $@ 2> $@


# Generate all plots and tables that summarize the results of the analyses. 
R/plots_tables.Rout: R/neut_AZ_full_M1_estwts.Rout R/spike_AZ_full_M1_estwts.Rout \
	R/neut_AZ_truncated_M2_estwts.Rout R/spike_AZ_truncated_M2_estwts.Rout
	Rscript --verbose R/plots_tables.R AZ > $@ 2> $@
