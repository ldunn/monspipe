import sys
import time
import ephem
import numpy as np
import matplotlib.pyplot as plt

def transit_correction(UTC_start, PSR_ra, PSR_dec, integration_time):

    # Geo-data for Molonglo
    # https://sites.google.com/a/utmostproject.org/internal/telescope-specifics/location
    # -35° 22' 14.5452" (-35.370707 deg), 149° 25' 28.7682" (149.424658 deg).
    utmost = ephem.Observer()
    utmost.lon = '149.424658' # deg
    utmost.lat = '-35.370707' # deg
    utmost.elevation = 732  # meters, obtained from google maps

    # relative strength of a source at the transit point
    amp = 1.0 

    # compute sigma for the NS-arm beam, for modeling as a Gaussian
    FWHP = 2.0  # width of NS-arm beam in the EW direction
    sigma = FWHP/(2.0*np.sqrt(2.0*np.log(2.0)))

    # remove third "-" from UTC string to get in right format for the ephem module
    UTC_start = UTC_start.split('-')
    UTC_start = UTC_start[0]+"-"+UTC_start[1]+"-"+UTC_start[2]+" "+UTC_start[3]

    #print("UTC_start = ", UTC_start)
    
    utmost.date = UTC_start
    utmost_lst = str(utmost.sidereal_time())
    utmost_lst = utmost_lst.split(":")
    utmost_lst = float(utmost_lst[0]) + float(utmost_lst[1])/60.0 + float(utmost_lst[2])/3600.0
    utmost_lst = np.around(utmost_lst,3)

    #print("Start LST = ",utmost_lst)

    #PSR_ra = PSR_ra.split(":")
    #PSR_ra = float(PSR_ra[0]) + float(PSR_ra[1])/60.0 + float(PSR_ra[2])/3600.0
    #PSR_ra = np.around(PSR_ra,3)

    #print("Source RA = ",PSR_ra,"hours")

    # meridian angle of the start of the observation (in degrees)
    theta_start = -(PSR_ra - utmost_lst)*15.0
    #print(theta_start)
    #print(PSR_dec)

    # offset to transit point in degrees - determined experimentally from data
    # this is the difference between the EM transit point and the local meridian (HA=0)
    theta_start -= 0.45

    # compute pulsar declination to nearest degree
    # this affects to transit rate through the EW primary beam as 1.0/cos(DEC)
    # meridian angle at the start of the observation in degrees 
    theta_start = theta_start * np.cos(PSR_dec*np.pi/180.0)

    # now use integration time in seconds to figure out meridian angle at end of obs (in degrees)
    theta_end = theta_start + (integration_time/3600.0)*15.0*np.cos(PSR_dec*np.pi/180.0)

    # set up an array of angles between theta_start and theta_end
    theta = np.linspace(theta_start,theta_end,101)

    # relative strength of a transiting and a tracked source
    # the primary beam response is modelled as a Gaussian for the transiting
    # source, and is constant for a tracked source
    flux = amp*np.exp(-(theta**2)/(2.0*sigma**2))
    const = amp*np.ones(len(theta))

    # figure out accumulated SNR of the transiting source over the start to end angles
    sum_of_squares = 0.0
    sn_cumulative_transit = np.zeros(len(flux))
    for i in range(len(flux)):
        sum_of_squares += flux[i]**2
        sn_cumulative_transit[i] += sum_of_squares
    sn_cumulative_transit = np.sqrt(sn_cumulative_transit)

    # and the accumulated SNR of the tracked source
    sum_of_squares = 0.0
    sn_cumulative_tracking = np.zeros(len(flux))
    for i in range(len(flux)):
        sum_of_squares += const[i]**2
        sn_cumulative_tracking[i] += sum_of_squares
    sn_cumulative_tracking = np.sqrt(sn_cumulative_tracking)

    # final SNR of the tracked source
    norm = sn_cumulative_tracking[-1]

    sn_cumulative_transit /= norm
    sn_cumulative_tracking /= norm

    #print("SNR of transiting source = ",sn_cumulative_transit[-1])
    #print("SNR of tracked source = ",sn_cumulative_tracking[-1])
    
    return sn_cumulative_transit[-1]/sn_cumulative_tracking[-1], PSR_ra - utmost_lst 

if __name__ == '__main__':

    # data for a test observation
    psr = sys.argv[1]
    UTC_start = sys.argv[2]
    integration_time = float(sys.argv[3])
    # seconds

    psr_ra = float(psr[1:3]) + float(psr[3:5])/60 # in hours
    if "-" in psr:
        dec = "-" + psr.split("-")[-1]
    elif "+" in psr:
        dec = "+" + psr.split("+")[-1]
    if len(dec) == 3:
        psr_dec = float(dec) # in degrees
    elif len(dec) == 5:
        psr_dec = float(f"{dec[:3]}")



    norm_factor, lst_offset = transit_correction(UTC_start, psr_ra, psr_dec, integration_time)

    print(norm_factor)
    print(lst_offset)
