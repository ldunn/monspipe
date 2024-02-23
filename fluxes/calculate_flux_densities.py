#!/usr/bin/env python

import psrchive
import subprocess
from transit_correction import transit_correction
import matplotlib.pyplot as plt
import numpy as np
import psycopg2
import sys
import json
def get_mean(x,y,N):
    xm = []
    ym = []
    for i in range(len(x)-N):
        xm = np.append(xm,np.mean(x[i:i+N]))
        ym = np.append(ym,np.mean(y[i:i+N]))
    return xm, ym

con = psycopg2.connect("dbname=utmost2d host=tooarrana2 user=ldunn")
cur = con.cursor()

calibrator_psr_list = "/fred/oz002/ldunn/UTMOST2D/code/calibrator_psrs.txt"

flux_correction_factor = 0.9367336183262432 # Calculated from flux_densities/flux_densities.ipynb, calibrated against Fabian's S843 values

with open(calibrator_psr_list, 'r') as f:
    calibrator_psrs = f.read().splitlines()

profile_dir = "/fred/oz002/ldunn/UTMOST2D/utmost-2d-profiles"

performance_min_cut = 0.35
performance_max_cut = 1.25

mjd_2021 = 59215
mjd_2022 = 59580
calibration_time = [mjd_2022 + 18.25, mjd_2022 + 73]
print(calibration_time)
calib_psrs_mean_snrs = {}

for psr in calibrator_psrs:

    if "-" in psr:
        dec = "-" + psr.split("-")[-1]
    elif "+" in psr:
        dec = "+" + psr.split("+")[-1]

    psr_ra = float(psr[1:3]) + float(psr[3:5])/60 # in hours

    if len(dec) == 3:
        psr_dec = float(dec) # in degrees
    elif len(dec) == 5:
        psr_dec = float(f"{dec[:3]}")


    cur.execute(f"SELECT snr,tim,date FROM observations WHERE name='{psr}' AND mjd > {calibration_time[0]} AND mjd < {calibration_time[1]} AND manual_exclude=false")
    calib_snrs = []
    for snr, tim, date in cur.fetchall():
        if snr < 7:
            continue
        tim_split = tim.split(' ')
        tobs = float(tim_split[tim_split.index('-tobs')+1])
        norm_snr = snr*np.sqrt(300/tobs)
        transit_corr, lst_offset = transit_correction(date, psr_ra, psr_dec, tobs)
        if np.abs(lst_offset) > 0.5:
            print("Observation is more than 0.5 hours away from PSR transit, skipping")
            continue
        calib_snrs.append(norm_snr/transit_corr)

    calib_psrs_mean_snrs[psr] = np.mean(calib_snrs)

print(calib_psrs_mean_snrs)
calib_psrs_running_mean = []

calib_psr_query = ""
for psr in calibrator_psrs[:-1]:
    calib_psr_query += f"name ='{psr}' OR "
calib_psr_query += f"name='{calibrator_psrs[-1]}'"

cur.execute(f"SELECT name,snr,tim,mjd FROM observations WHERE ({calib_psr_query}) ORDER BY mjd")

calib_psrs_snrs = []

for psr, snr, tim, mjd in cur.fetchall():
    if snr < 7: 
        continue

    tim_split = tim.split(' ')
    tobs = float(tim_split[tim_split.index('-tobs')+1])
    norm_snr = snr*np.sqrt(300/tobs)

    calib_snr = norm_snr/calib_psrs_mean_snrs[psr]
    calib_psrs_snrs.append((mjd, calib_snr))

calib_psrs_snrs = np.array(calib_psrs_snrs)
print(calib_psrs_snrs)
performance_mask = np.logical_and(calib_psrs_snrs[:,1] > performance_min_cut, calib_psrs_snrs[:,1] < performance_max_cut)

calib_psrs_snrs = calib_psrs_snrs[performance_mask, :]

mean_times, mean_snrs = get_mean([x[0] for x in calib_psrs_snrs], [x[1] for x in calib_psrs_snrs], 15)

print(len(mean_snrs))

psr = sys.argv[1]
date = sys.argv[2]

Tsky_file = "/fred/oz002/ldunn/UTMOST2D/code/psr_tskys.dat"
Tskys = np.genfromtxt(Tsky_file, dtype=None)
Tskys = dict([(entry[0].decode('UTF-8'), entry[3]) for entry in Tskys])

if psr not in Tskys:
    print(f"Could not find a Tsky for {psr}. Aborting.")
    exit(1)

Tsky = Tskys[psr]
print(f"Tsky: {Tsky}")

Tsys = 170
gain = 0.227
SEFD = (Tsys + Tsky)/gain
#SEFD = 750
MOLONGLO_LAT = -35.4

if "-" in psr:
    dec = "-" + psr.split("-")[-1]
elif "+" in psr:
    dec = "+" + psr.split("+")[-1]

psr_ra = float(psr[1:3]) + float(psr[3:5])/60 # in hours

if len(dec) == 3:
    psr_dec = float(dec) # in degrees
elif len(dec) == 5:
    psr_dec = float(f"{dec[:3]}")

print(f"Pulsar declination: {psr_dec}")

meridian_correction = 1.0/np.cos((psr_dec - MOLONGLO_LAT) * np.pi/180)

print(f"SELECT name,snr,tim,date,mjd FROM observations WHERE name='{psr}' AND date='{date}'")
cur.execute(f"SELECT name,snr,tim,date,mjd FROM observations WHERE name='{psr}' AND date='{date}'")
for psr, snr, tim, date, mjd in cur.fetchall():
    if snr < 7:
        print("SNR < 7. Not calculating a flux density.")
        flux = None
    else:
        tim_split = tim.split(' ')
        tobs = float(tim_split[tim_split.index('-tobs')+1])
        norm_snr = snr*np.sqrt(300/tobs)

        nearest_time_idx = np.argmin([np.abs(t-mjd) for t in mean_times])
        
        transit_corr, lst_offset = transit_correction(date, psr_ra, psr_dec, tobs)

        if np.abs(lst_offset) > 0.5:
            print("Observation is more than 0.5 hours away from PSR transit, skipping")
            continue
        print(f"Transit correction: {transit_corr}")

        corrected_snr = norm_snr/mean_snrs[nearest_time_idx]/transit_corr
        profile_archive = f"{profile_dir}/{psr}.std"
        psrstat_cmd = f"psrstat -c width=transitions:{{fmax=0.5}} -c width {profile_archive}"
        psrstat_out = subprocess.check_output(psrstat_cmd, encoding="UTF-8", shell=True)
        width = float(psrstat_out.split("=")[-1])

        flux = SEFD*np.sqrt(width/(1-width))*corrected_snr/np.sqrt(300*2*45e6) # BW is 50MHz - 2.5MHz auto-zapped on each edge
        flux *= meridian_correction
        flux /= flux_correction_factor
        print(flux)

    results = {
            "percent_rfi_zapped": None,
            "dm": None,
            "dm_err": None,
            "dm_epoch": None,
            "dm_chi2r": None,
            "dm_tres": None,
            "rm": None,
            "rm_err": None,
            "sn": snr,
            "flux": flux,
        }
    with open(f"/fred/oz002/ldunn/meertime_dataportal/data/post/{psr}/{date}/results.json", "w") as f:
        json.dump(results, f, indent=1)
