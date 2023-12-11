This repository contains the scripts used to process the observations from the UTMOST-NS pulsar timing programme.
The code for the HMM glitch detection pipeline can be found [here](https://github.com/ldunn/glitch_hmm).

A very brief description of the files included:

* `process_single_obs_nt.sh` - this shell script takes the raw archive files produced at Molonglo and processes them into RFI-cleaned archives using `xprof`, as well as producing a ToA
* `make_plots.sh` - this is called by `process_single_obs_nt.sh` to produce some useful plots from both the raw and cleaned archives using PSRCHIVE
* `prepare_ingest.py` - this script runs the HMM glitch detection step, and produces a plot of the timing residuals via `libstempo`
* `phasing.log` - this file records all of the different phasing configurations of the instrument - this is how we break up observing "sessions" for the UTMOST-NS data
* `hmm_default.ini` - this is the configuration file for the HMM glitch detection step
* `fluxes/calculate_flux_densities.py` - this script calculates an estimate of the pulsar flux density for a given observation. It currently relies on access to an internal database, but all of the information used is contained in the archive files.
* `fluxes/transit_correction.py` - this contains a helper function to calculate the transit correction factor to be applied to each observation when calculating the flux density
* `fluxes/calibrator_psrs.txt` - this file contains a list of the pulsars used to calibrate the telescope sensitivity curve and hence the flux density measurements
