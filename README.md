# Cottrell2026
Repository of code used for Cottrell et al 2026 (PNAS) "HES1 oscillations are required for cell cycle reentry in oestrogen receptor–positive breast cancer cells".

Time-series were pre-processed in R using time-series_preprocessing.Rmd and periodicity was estimated in MatLab using Lomb-scargle periodogram (LSP_code.m) and Autocorrelation Function (Autocorrelation_code.m). Time-series were further analysed by identifying peaks and dips using Peaks_dips_detection.m, this allowed calculation of amplitudes, fold-changes and relative timing of dynamioc features from co-imaging experiments.
