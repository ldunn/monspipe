#!/usr/bin/env python

import glob
import logging
import time
import os
import json
import getpass
from datetime import datetime, timedelta, timezone
import tempfile
import astropy
import psrchive
import libstempo
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import configparser
import pathlib
import subprocess
import pulsar_hmm.HMM as HMM

from tables import *
from graphql_client import GraphQLClient
from util import header, ephemeris
from util import time as util_time

RESULTS_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/post'
FOLDING_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/post'
IMAGES_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/images'
PROFILES_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/profiles'
PARS_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/pars'
EXCLUDED_TOAS = '/fred/oz002/ldunn/meertime_dataportal/data/excluded_toas.dat'
CALIBRATIONS_FILE = '/fred/oz002/ldunn/meertime_dataportal/pipeline/phasing.log'
HMM_CONFIG = '/fred/oz002/ldunn/meertime_dataportal/pipeline/hmm_default.ini'
HMM_WORKING_DIR = '/fred/oz002/ldunn/meertime_dataportal/data/hmm_working_dir'
FLUX_DENSITY_SCRIPT = '/fred/oz002/ldunn/UTMOST2D/code/calculate_flux_densities.py'
DT_FMT = "%Y-%m-%d-%H:%M:%S"
DT_FMT_GRAPHQL = "%Y-%m-%dT%H:%M:%S%z"

SNR_THRESH=7

images_dict = {"fy": "phase-time",
               "gt": "phase-freq",
               "dftp": "profile-int",
              }

def get_calibration(utc):
    calibrations = "None"
    with open(CALIBRATIONS_FILE) as f:
        for line in f.readlines():
            if line[0] == '#':
                continue
            
            date = utc.strptime(line.split(' ')[0], '%Y-%m-%d').replace(tzinfo=timezone.utc)
            delta = date - utc
            if delta.total_seconds() > 0:
                break

            calibration = line.split(' ')[1].strip()

    return calibration

def parse_toa(toa_text):
    split = toa_text.split(' ')
    freq = float(split[1])
    mjd = split[2]
    uncertainty = float(split[3])
    site = split[4][0]
    flags = {}
    for flag, val in zip(split[5:-1:2], split[6::2]):
        flags[flag[1:]] = val

    return json.dumps(flags), freq, mjd, site, float(flags['snr']), uncertainty

def get_timfile(source, utc_start_dt, utc_end_dt, snr_thresh=SNR_THRESH, excluded_toas=None, utc_keep_dt=None):
    all_obs = glob.glob(f"{RESULTS_DIR}/{source}/2*")
    tims = []

    for obs in all_obs:
        
        utc = obs.split("/")[-1]
        utc_dt = datetime.strptime(utc, DT_FMT).replace(tzinfo=timezone.utc)
        if excluded_toas is not None:
            print(utc)
            print(utc in excluded_toas[source])
            if utc in excluded_toas[source] and utc_dt != utc_keep_dt:
                continue

        if not (utc_dt <= utc_end_dt and utc_dt >= utc_start_dt):
            continue
        
        toa_file = f"{obs}/toa.tim"
        if not os.path.isfile(toa_file):
            continue

        tim_line = open(toa_file).read().strip()
        _, _, _, _, snr, _ = parse_toa(tim_line)
        if snr < snr_thresh and utc_dt != utc_keep_dt:
            continue

        tims.append(tim_line)

    f = tempfile.NamedTemporaryFile(mode='w')
    print("MODE 1\nFORMAT 1", file=f)
    for tim in tims:
        print(tim, file=f)
    timfile = f.name
    f.flush()
    return f, tims


def setup_hmm(par, tim, config):
    psr = libstempo.tempopulsar(parfile=par, timfile=tim)
    print(par)
    print(tim)
    print(psr['F0'].val)
    print(psr.toaerrs)

    freq_min = float(config['doi']['freq_min'])
    freq_max = float(config['doi']['freq_max'])
    dfreq = float(config['doi']['dfreq'])
    fdot_min = float(config['doi']['fdot_min'])
    fdot_max = float(config['doi']['fdot_max'])
    dfdot = float(config['doi']['dfdot'])

    freqs = np.arange(freq_min, freq_max, dfreq)
    fdots = np.arange(fdot_min, fdot_max, dfdot)
    if len(fdots) % 2 == 0:
        fdots = np.linspace(fdot_min, fdot_max, len(fdots)+1)
        dfdot = np.diff(fdots)[0]
        print(f'Odd number of fdot bins required, adjusting dfdot. New dfdot = {dfdot}')

    min_toa_gap = 0
    mjd_range = None
    psr_name = None
    red_idx = None
    red_amp = None
    red_comp = None
    efac = 1
    equad = -100

    with open(par,'r') as f:
        for line in f.readlines():
            split = line.split()
            if split[0].strip() == 'PSRJ':
                psr_name = split[-1].strip()
                print(psr_name)
            elif split[0].strip() == 'TNRedAmp':
                red_amp = 10**float(split[-1].strip())
            elif split[0].strip() == 'TNRedGam':
                red_idx = float(split[-1].strip())
            elif split[0].strip() == 'TNRedC':
                red_comp = float(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobalef':
                efac = float(split[-1].strip())
            elif split[0].strip().lower() == 'tnglobaleq':
                equad = float(split[-1].strip())

    tn_params = {'red_amp': red_amp, 'red_idx': red_idx, 'red_comp': red_comp, 'efac': efac, 'equad': equad}

    if 'tn' in config:
        if 'sigma' in config['tn']:
            print('Overriding default value of sigma')
            sigma = float(config['tn']['sigma'])
    else:
        print(psr)
        mean_z = np.mean(np.diff(sorted(psr.toas())))*86400
        sigma = max(dfdot/np.sqrt(mean_z), 1e-21)
    hmm = HMM.HMM.from_tempo2(par, tim, freqs, fdots, None, [], min_toa_gap = min_toa_gap, mjd_range=mjd_range)
    return hmm, psr, sigma, tn_params

def save_hmm_files(hmm, working_prefix):
    np.savetxt(f"{working_prefix}kappas.dat", hmm.kappas)
    np.savetxt(f"{working_prefix}freqs.dat", hmm.freqs)
    np.savetxt(f"{working_prefix}fdots.dat", hmm.fdots)
    np.savetxt(f"{working_prefix}zs.dat", hmm.zs)

def do_hmm_psr(hmm, sigma, config, out_prefix, working_prefix, extra_matlab_cmd=None):
    matlab_wrapper = config['matlab']['matlab_wrapper']
    matlab_path = config['matlab']['matlab_path']
    #out_prefix = config['out']['out_prefix']
    #working_prefix = config['general']['working_prefix']

    matlab_cmd = "global f_fiducial fd_fiducial fdd_fiducial sigma kappa_per_toa;"
    matlab_cmd += f"f_fiducial = {hmm.f_fiducial};"
    matlab_cmd += f"fd_fiducial = {-hmm.fd_fiducial};"
    matlab_cmd += f"fdd_fiducial = {hmm.fdd_fiducial};"
    matlab_cmd += f"sigma = {sigma};"
    matlab_cmd += f"out_prefix = '{out_prefix}';"
    matlab_cmd += f"working_prefix = '{working_prefix}';"
    matlab_cmd += f"matlab_path = '{matlab_path}';"
    if extra_matlab_cmd:
        matlab_cmd += extra_matlab_cmd
    matlab_cmd += f"run('{matlab_wrapper}');"
    print(matlab_cmd)
    cmd = f"matlab -nosplash -nodesktop -r \"{matlab_cmd} exit\"" 
    p = subprocess.run(cmd, shell=True)
    print(f"MATLAB process returned error code {p.returncode}")
    return p.returncode == 0

# Rewrite .tim for HMM to contain PETs (ToAs in pulsar reference frame)
def save_hmm_tim(toas, par, working_prefix):
    tim = "FORMAT 1\nMODE 1"
    for toa in toas:
        tim += "\n"+toa

    bats_fname = working_prefix + ".bats"
    with open(bats_fname, "w") as f:
        print(tim, file=f)

    tempo2_cmd = f"tempo2 -output general2 -f {par} {bats_fname} -nofit -s \
                \"fake {{freq}} {{pet}} {{err}} BAT {{flags}} -phaseJ {{phaseJ}}\\n\""

    pets_output = subprocess.check_output(tempo2_cmd, shell=True).decode('utf-8')
    tim_lines = [line for line in pets_output.split('\n') if 'phaseJ' in line]
    tim = "FORMAT 1\nMODE 1"
    for toa in tim_lines:
        tim += "\n" + toa

    tim_fname = working_prefix + ".tim"
    with open(tim_fname, "w") as f:
        print(tim, file=f)
        f.flush()

    return working_prefix + ".tim"

def run_hmm(source, utc_start_dt, utc_end_dt, config, excluded_toas=None):
    
    timfile, tims = get_timfile(source, utc_start_dt, utc_end_dt, excluded_toas=excluded_toas)
    if len(tims) < 5:
        logging.warning(f"Number of available ToAs between {utc_start_dt.strftime(DT_FMT)} and {utc_end_dt.strftime(DT_FMT)} for {source} is less than 5. Aborting HMM run.")
        return {}

    mjds = [tim.split(" ")[2] for tim in tims]

    parfile = f"{PARS_DIR}/{source}.par"
    if not os.path.isfile(parfile):
        logging.warn(f"Do not have par file for {psr}! Aborting HMM run.")
        return {}

    if len(mjds) > 30:
        timfile.close()
        sort = np.argsort(mjds)
        print(sort[-30:])
        tims = np.array(tims)[sort[-30:]]
        print(tims)
#        timfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
#        timfile.flush()

    hmm_tim = save_hmm_tim(tims, parfile, timfile.name)
    
    print(hmm_tim)
    out_prefix = f"{RESULTS_DIR}/{source}/hmm_runs/{utc_end_dt.strftime(DT_FMT)}/{source}_"
    plot_prefix = f"{IMAGES_DIR}/{source}/hmm_runs/{utc_end_dt.strftime(DT_FMT)}/{source}_"
    working_prefix = f"{HMM_WORKING_DIR}/{source}_{utc_end_dt.strftime(DT_FMT)}/{source}_"

    pathlib.Path(out_prefix).parent.absolute().mkdir(exist_ok=True, parents=True)
    pathlib.Path(plot_prefix).parent.absolute().mkdir(exist_ok=True, parents=True)
    pathlib.Path(working_prefix).parent.absolute().mkdir(exist_ok=True, parents=True)

    
    hmm, _, sigma, tn_params = setup_hmm(parfile, hmm_tim, config)

    save_hmm_files(hmm, working_prefix)
    if not do_hmm_psr(hmm, sigma, config, out_prefix, working_prefix):
        return {}

    make_hmm_plots(hmm, out_prefix, plot_prefix)
    hmm_results = {}
    hmm_results["utc_start"] = utc_start_dt.strftime(DT_FMT)
    hmm_results["utc_end"] = utc_end_dt.strftime(DT_FMT)
    hmm_results["lnBF_thresh"] = 10**0.5
    hmm_results["lnBFs"] = list(np.atleast_2d(np.loadtxt(f"{out_prefix}bfs.dat"))[0,:])
    hmm_results["glitch_cand"] = any([x > hmm_results["lnBF_thresh"] for x in hmm_results["lnBFs"]])
    with open(f"{out_prefix}results.json", "w") as f:
        json.dump(hmm_results, f)
    return hmm_results
    
def make_hmm_plots(hmm, out_prefix, plot_prefix):
    freqs = hmm.freqs
    fdots = hmm.fdots
    zs = hmm.zs
    f_posterior = np.loadtxt(f"{out_prefix}f_posterior.dat")
    plt.imshow(np.exp(f_posterior), aspect='auto')
    plt.title("Frequency posterior")
    plt.colorbar(label="Probability density")
    plt.xlabel('ToA index')
    plt.ylabel('Frequency bin')
    plt.gca().invert_yaxis()
    plt.savefig(f"{plot_prefix}f_posterior.png")
    plt.clf()

    f_path = np.loadtxt(f"{out_prefix}f_path.dat", dtype=np.float)
    plt.plot(np.cumsum(zs)/86400, f_path)
    plt.title("Maximum a-posteriori frequency path")
    plt.ylim(min(freqs), max(freqs))
    plt.xlabel('Days since first ToA')
    plt.ylabel('Frequency (Hz)')
    plt.savefig(f"{plot_prefix}f_path.png")
    plt.clf()

    fdot_posterior = np.loadtxt(f"{out_prefix}fdot_posterior.dat")
    plt.imshow(np.exp(fdot_posterior), aspect='auto')
    plt.title("Frequency derivative posterior")
    plt.colorbar(label="Probability density")
    plt.xlabel('ToA index')
    plt.ylabel('Frequency derivative bin')
    plt.gca().invert_yaxis()
    plt.savefig(f"{plot_prefix}fdot_posterior.png")
    plt.clf()


    fdot_path = np.loadtxt(f"{out_prefix}fdot_path.dat", dtype=np.float)
    print(fdot_path)
    plt.plot(np.cumsum(zs)/86400, fdot_path)
    plt.title("Maximum a-posteriori frequency derivative path")
    plt.ylim(min(fdots), max(fdots))
    plt.xlabel('Days since first ToA')
    plt.ylabel('Frequency derivative (Hz/s)')
    plt.savefig(f"{plot_prefix}fdot_path.png")
    plt.clf()

    bfs = np.atleast_2d(np.loadtxt(f"{out_prefix}bfs.dat"))
    for i in range(0, bfs.shape[0]):
        plt.plot(np.cumsum(zs)/86400, bfs[i, :])
        plt.title("ln Bayes factor glitch vs. no glitch")
        plt.xlabel('Days since first ToA')
        plt.ylabel('ln Bayes factor')
        plt.savefig(f"{plot_prefix}bfs_{i}.png")
        plt.clf()

    return


def make_residual_plot(source, obs_dt, utc_end_dt, excluded_toas=None, highlight_utc_dt=None):
    
    utc_start_dt = datetime.utcfromtimestamp(0).replace(tzinfo=timezone.utc)
    timfile_handle, _ = get_timfile(source, utc_start_dt, utc_end_dt, excluded_toas=excluded_toas, utc_keep_dt=highlight_utc_dt)
    timfile = timfile_handle.name
    psr = libstempo.tempopulsar(parfile=f"{PARS_DIR}/{source}.par", timfile=timfile)
    matplotlib.use('Agg')
    toas = sorted(psr.toas())
    residuals = np.asarray([x for _, x in sorted(zip(psr.toas(), psr.residuals()))])
    toaerrs = np.asarray([x for _, x in sorted(zip(psr.toas(), psr.toaerrs))])
    
    plt.errorbar(toas, residuals*psr['F0'].val, toaerrs/1e6*psr['F0'].val, fmt='.', capsize=3, alpha=0.2)
    plt.ylabel('res [turns]')
    plt.xlabel("MJD")
    phase_ax = plt.gca()
    phase_ax.set_yticklabels([f"{x:.2g}" for x in phase_ax.get_yticks().tolist()])
    time_ax = phase_ax.twinx()
    plt.errorbar(toas, residuals*psr['F0'].val, toaerrs/1e6*psr['F0'].val, fmt='.', capsize=3, alpha = 0.2)

    if highlight_utc_dt is not None:
        highlight_mjd = astropy.time.Time(highlight_utc_dt, scale='utc').mjd
        highlight_idx = (np.abs(toas - highlight_mjd)).argmin()
        plt.errorbar(toas[highlight_idx], residuals[highlight_idx]*psr['F0'].val, toaerrs[highlight_idx]/1e6*psr['F0'].val, fmt='*', capsize=3)

    time_ax.set_ylabel('Residuals [sec]')
    phase_ax.set_ylabel('Residuals [turns]')
    time_axlabels_floats = [float(x.get_text())/psr['F0'].val for x in phase_ax.get_yticklabels()]
    time_ax.set_yticklabels([f"{label:.2E}" for label in time_axlabels_floats])
    rms_res = np.sqrt(np.mean(psr.residuals()**2))
    plt.title(f"rms residuals: {rms_res:.3e} [s]")
    plt.tight_layout()
    
    plt.savefig(f"{IMAGES_DIR}/{source}/{obs_dt.strftime(DT_FMT)}/{source}_{obs_dt.strftime(DT_FMT)}_tempo2.png")

    timfile_handle.close()

def parse_excluded_toas():
    excluded_toas = {}
    for line in open(EXCLUDED_TOAS, "r").readlines():
        split = line.strip().split(",")
        if len(split) > 1:
            excluded_toas[split[0]] = split[1:]

    return excluded_toas

def main(source, utc_start):

    excluded_toas = parse_excluded_toas()
    results_dir = f"{RESULTS_DIR}/{source}/{utc_start}"
    
    obs_header = f"{results_dir}/obs.header"
    if not os.path.isfile(obs_header):
        error = Exception(f"{obs_header} not found.")
        logging.error(str(error))
        raise error

    try:
        hdr = header.Header(obs_header)
    except Exception as error:
        logging.error(str(error))
        raise error

    hdr.set("TELESCOPE", "MONS")
    hdr.parse()

    archive_raw_FT_path = f"{results_dir}/{source}_{utc_start}.raw.FT"
    if not os.path.isfile(archive_raw_FT_path):
        error = Exception(f"{archive_raw_FT_path} not found.")
        logging.error(str(error))
        raise error

    archive_raw_FT = psrchive.Archive_load(archive_raw_FT_path)

    archive_clean_FT_path = f"{results_dir}/{source}_{utc_start}.FT"
    if not os.path.isfile(archive_clean_FT_path):
        error = Exception(f"{archive_clean_FT_path} not found.")
        logging.error(str(error))
        raise error

    archive_clean_FT = psrchive.Archive_load(archive_clean_FT_path)


    literal = False
    quiet = True

    utc_start_dt = datetime.strptime(utc_start, "%Y-%m-%d-%H:%M:%S").replace(tzinfo=timezone.utc)
    
    duration = archive_raw_FT.get_Integration(0).get_duration()
    suspect = False
    comment = ""

    raw_eph = ephemeris.Ephemeris()
    raw_eph.load_from_archive_as_str(archive_raw_FT_path)

    clean_eph = ephemeris.Ephemeris()
    clean_eph.load_from_archive_as_str(archive_clean_FT_path)

    created_at = util_time.get_current_time()
    created_by = getpass.getuser()
    valid_from = util_time.get_time(0)
    valid_to = util_time.get_time(4294967295)
    comment = "Created by tables.ephemeris.new"

    location = f"{RESULTS_DIR}/{source}/{utc_start}/"

    raw_results = json.dumps({"snr": archive_raw_FT.get_Integration(0).get_Profile(0,0).snr()})

    location = f"{RESULTS_DIR}/{source}/{utc_start}/"
    clean_flux_cmd = [FLUX_DENSITY_SCRIPT, source, utc_start]
    print(' '.join(clean_flux_cmd))
    toa_path = f"{RESULTS_DIR}/{source}/{utc_start}/toa.tim"
    snr = None
    if os.path.isfile(toa_path):
        with open(toa_path) as f:
            toa_text = f.read().strip()
        
        flags, freq, mjd, site, snr, uncertainty = parse_toa(toa_text)

    if utc_start not in excluded_toas[source]:
        try:
            flux = float(subprocess.check_output(clean_flux_cmd).decode("UTF-8").strip().split("\n")[-1].strip())
            clean_results = json.dumps({"snr": snr if snr is not None else archive_clean_FT.get_Integration(0).get_Profile(0,0).snr(), "flux": flux*1000, "zap_frac": 0})
            print(clean_results)
        except:
            print(f"Flux density calculation failed.")
            clean_results = json.dumps({"snr": snr, "zap_frac": 0})
    else:
        clean_results = json.dumps({"snr": snr, "zap_frac": 0})


    template_archive_path = f"{PROFILES_DIR}/{source}.std"
    template_archive = psrchive.Archive_load(template_archive_path)
    freq = template_archive.get_centre_frequency()
    bw = template_archive.get_bandwidth()
    
    toa_path = f"{RESULTS_DIR}/{source}/{utc_start}/toa.tim"
    if os.path.isfile(toa_path):
        with open(toa_path) as f:
            toa_text = f.read().strip()
        
        flags, freq, mjd, site, snr, uncertainty = parse_toa(toa_text)

        comment = {"tim": toa_text}
        print(utc_start)
        print(excluded_toas[source])
        print(utc_start in excluded_toas[source])
    else:
        logging.info("Did not find a ToA for this observation")
   
    if snr > SNR_THRESH and utc_start not in excluded_toas[source]:
        hmm_utc_end_dt = utc_start_dt
        hmm_utc_start_dt = hmm_utc_end_dt + timedelta(days=-10000)
        hmm_config = configparser.ConfigParser()
        hmm_config.read(HMM_CONFIG)
        config_to_dump = {}
        for section in hmm_config.sections():
            for key, val in hmm_config.items(section):
                config_to_dump[key] = val
        hmm_results = run_hmm(source, hmm_utc_start_dt, hmm_utc_end_dt, hmm_config, excluded_toas=excluded_toas)
    else:
        hmm_processing_id = None

    make_residual_plot(source, utc_start_dt, datetime.strptime("2024-01-01-00:00:00", DT_FMT).replace(tzinfo=timezone.utc), excluded_toas=excluded_toas,highlight_utc_dt=utc_start_dt)
    
if __name__ == "__main__":
    import argparse

    import faulthandler
    faulthandler.enable()
    parser = argparse.ArgumentParser(description="Ingest PTUSE fold mode observation")
    parser.add_argument("source", type=str, help="source of the obs")
    parser.add_argument("utc_start", type=str, help="utc_start of the obs")
    args = parser.parse_args()
    #logging.basicConfig(level=logging.DEBUG)
    main(args.source,args.utc_start)
