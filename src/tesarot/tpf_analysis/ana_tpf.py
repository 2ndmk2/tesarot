import numpy as np
import lightkurve as lk
from tesarot.utils import utils
import os
from astropy.table import Table
from tesarot.lc_analysis import lc_ana


def get_tess_tpf(target, mission = "TESS", exptime = 120):
    search_result = lk.search_targetpixelfile(target, mission=mission, exptime = exptime)
    tpfs = search_result.download_all()
    return tpfs, search_result

def save_search_result_tpf(search_result, target, outdir):
    table = search_result.table
    file_name =os.path.join(outdir, "search_resut_tpf_%s.csv" % target)
    table.write(file_name, format ="csv")

def save_tpf(tpfs, search_result,  outdir):

    file_heads = utils.make_output_name_for_tpfs(search_result)
    for (i, tpf) in enumerate(tpfs):
        file_name = os.path.join(outdir, file_heads [i]) + ".fits"
        tpf.to_fits(file_name )
import lightkurve as lk

def load_tpf(search_result, outdir):

    file_heads = utils.make_output_name_for_tpfs(search_result)
    tpfs = []
    for (i, file_name) in enumerate(file_heads):
        file_name = os.path.join(outdir, file_heads [i]).replace(" ","")+ ".fits"
        tpfs.append(lk.read(file_name))
    return tpfs

def load_data_tpf(outdir, target):
    result_file_name =os.path.join(outdir, "search_resut_tpf_%s.csv" % target)
    table_result = Table.read(result_file_name, format='csv', fast_reader=False)
    search_result = lk.SearchResult(table = table_result )
    tpfs = load_tpf(search_result, outdir)
    return tpfs, search_result
    
def compute_centroid(image, x_cen, y_cen, wd):
    center_x = 0
    center_y = 0
    flux_all = 0
    i_min = int(y_cen - wd + 1)
    i_max = int(y_cen + wd + 1)
    j_min = int(x_cen - wd+ 1)
    j_max = int(x_cen + wd+ 1)    

    for i in range(i_min, i_max):
        for j in range(j_min, j_max):
            center_x += j* image[i][j]
            center_y += i* image[i][j]
            flux_all += image[i][j]
    return center_x/flux_all, center_y/flux_all

def make_difference_image_for_tpfs(tpfs, period_inputs = None, aperture_mask= "default", filter_length = 720):
    diff_images = []
    folded_lcs = []
    folded_bin_lcs = []
    avg_images = []
    periods = []

    for (i, tpf) in enumerate(tpfs):

        print ("diff_image:%d" % i)
        lc = tpf.to_lightcurve(aperture_mask=aperture_mask)
        lc_reduced = lc_ana.reduce_lc_for_periodogram(lc, filter_length)
        pg = lc_reduced.to_periodogram(oversample_factor=1)

        if period_inputs is None:
            period = pg.period_at_max_power
        else:
            period = period_inputs[i]
        folded = lc_reduced.fold(period, epoch_time=0.5)
        folded_bin = folded.bin(time_bin_size=0.01)
        full_phase_range = folded_bin.phase[-1].value - folded_bin.phase[0].value
        tolerance = 0.1 * full_phase_range
        min_phase = folded_bin.time[np.argmin(folded_bin.flux)].value
        max_phase = folded_bin.time[np.argmax(folded_bin.flux)].value
        min_timestamps = folded.time_original[np.where((folded.time > min_phase - tolerance)
                                                     & (folded.time < min_phase + tolerance))].value
        max_timestamps = folded.time_original[np.where((folded.time > max_phase - tolerance)
                                                     & (folded.time < max_phase + tolerance))].value  
        one_quarter_minima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in min_timestamps]
        one_quarter_maxima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in max_timestamps]    
        avg_image = np.nanmean(tpf.flux.value, axis=0)
        diff_image = np.abs(np.nanmean(one_quarter_maxima, axis=0) - np.nanmean(one_quarter_minima, axis=0))

        diff_images.append(diff_image)
        avg_images.append(avg_image)
        folded_lcs.append(folded)
        folded_bin_lcs.append(folded_bin)
        periods.append(period)

    return avg_images, diff_images, periods, folded_lcs, folded_bin_lcs

