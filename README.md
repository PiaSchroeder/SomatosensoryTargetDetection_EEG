# SomatosensoryTargetDetection_EEG
Matlab scripts to analyse data from an EEG study on somatosensory target detection.

Requires SPM12 software package: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

## Scripts
### Behavioural analysis:
* soma_eeg_check_pfs.m
* fit_logistic.m
* load_logs.m
* soma_eeg_behaviour.m

### Bayesian 1st level GLMs:
* soma_eeg_bayesglm_erp.m
* bayesglm_sensors.m

### Bayesian Model Selection:
* soma_eeg_bms_erp.m
* spm_BMS_family.m
* spm_dirichlet_exceedance_fam.m

### Plot results:
* soma_eeg_compute_int_erps.m
* soma_eeg_bms_results.m
* soma_eeg_compute_intmatched_det_erps.m
* soma_eeg_plot_topos.m
* soma_eeg_channel_count.m
