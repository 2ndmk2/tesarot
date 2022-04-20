import numpy as np
import lightkurve as lk
from tesarot.utils import utils
import os
from astropy.table import Table
from tesarot.lc_analysis import lc_ana

def make_timestamps_for_differential_imaging(lc, period):
    """ make timestamps for peaks & vallies of for periodic signals

    Params:
        lc: Lightcurve object
        period: Period

    Returns:
        min_timestamps: Arrays of time for vallies
        max_timestamps: Arrays of time for peaks
        folded: Folded lc
        folded_bin: Binned folex lc

    """


    folded = lc.fold(period, epoch_time=0.5)
    folded_bin = folded.bin(time_bin_size=0.01)
    full_phase_range = folded_bin.phase[-1].value - folded_bin.phase[0].value
    tolerance = 0.1 * full_phase_range
    min_phase = folded_bin.time[np.argmin(folded_bin.flux)].value
    max_phase = folded_bin.time[np.argmax(folded_bin.flux)].value
    min_timestamps = folded.time_original[np.where((folded.time > min_phase - tolerance)
                                                 & (folded.time < min_phase + tolerance))].value
    max_timestamps = folded.time_original[np.where((folded.time > max_phase - tolerance)
                                                 & (folded.time < max_phase + tolerance))].value  

    return min_timestamps, max_timestamps, folded, folded_bin

def make_diff_image(tpf, min_timestamps, max_timestamps, window_length=6001):
    """ Differential imaging based on tpf file

    Params:
        tpf: Target Pixel File object
        min_timestamps: Arrays of time for vallies
        max_timestamps: Arrays of time for peaks
        window_length (odd int): Length for filter in lk module 
    Returns:
        av_image: Averaged image
        diff_image: Differential image

    """

    nx, ny = np.shape(tpf.flux.value[0])
    diff_image = np.zeros((nx, ny))
    av_image = np.zeros((nx, ny))
    for k in range(nx * ny):
        x, y = np.unravel_index(k, (nx, ny))
        aperture_mask_now = np.zeros((nx, ny))
        aperture_mask_now[x,y] = 1
        lc = tpf.to_lightcurve(aperture_mask=aperture_mask_now )
        lc, trend_lc = lc.remove_outliers().flatten(window_length= window_length, polyorder=2, return_trend= True)
        lc.flux = lc.flux * np.median(trend_lc.flux)
        one_quarter_minima = [f for (f, t) in zip(lc.flux.value, lc.time.value) if t in min_timestamps]
        one_quarter_maxima = [f for (f, t) in zip(lc.flux.value, lc.time.value) if t in max_timestamps]  
        diff_flux = np.abs(np.nanmean(one_quarter_maxima) - np.nanmean(one_quarter_minima))
        diff_image[x, y]= diff_flux
        av_image[x,y] = np.nanmean(lc.flux.value)
    return av_image, diff_image


