EJT 2025

zoo_projections_train_RF: R code to train the NE Atlantic copepod zooplankton Random Forest (RF) machine learning models used to project zooplankton abundance under future climate change for Tyldesley et al 2025 Projected declines in zooplankton energy supporting Northeast Atlantic ecosystems [submitted]

scripts/pre_processing/clean_up_data.R: tidies zed_and_amm7.csv which contains Continuous Plankton Recorder samples matched in time and space with explanatory variables derived from the Copernicus NW European Shelf phys-bgc ocean reanalysis.

scripts/main/train_rf_final_dec2025.R: trains the RF models

scripts/main/predict_rf_projns.R: predicts the trained models on the bias corrected monthly projection ensemble derived from the RECICLE and HADGEM regional climate projections. Too large to host here but see MATLAB code repo for how to derive them: https://github.com/EmmaTyldesley/zoo_projns_paper_final_code.git

scripts/sensitivity_analysis: contains R scripts to optimise the RF parameters and check for ability to extrapolate in space and time

