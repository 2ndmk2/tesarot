import numpy as np
import lightkurve as lk
from tesarot.utils import utils
import os
from astropy.table import Table
from tesarot.lc_analysis import lc_ana
from tesarot.tpf_analysis import differential_imaging
from astropy.coordinates import SkyCoord


def get_tess_tpf(target, mission = "TESS", exptime = 120):

    """ Get target pixel files and search information
        Args:
            target: target name
            mission: mission name
            exptime: exposure time
        Returns:
            tpfs_arr: Arrays of target pixel files
            search_result: Table containing search result for data 

    """

    if exptime is None:
        search_result = lk.search_targetpixelfile(target, mission=mission)
    else:
        search_result = lk.search_targetpixelfile(target, mission=mission, exptime = exptime)
    tpfs = search_result.download_all()
    tpfs_arr = []

    if tpfs is None:
        print("no data:",target)
        return None, None

    for tpf in tpfs:
        tpfs_arr.append(tpf[tpf.quality<15])
    return tpfs_arr, search_result

def save_search_result_tpf(search_result, target, outdir):
    """ Save search result information for tpfs.

        Args:
            search_result: Table containing search result for data 
            target: Target name
            out_dir: Directory for output
        Returns:
            wcs_arr: Arrays of wcs information 

    """

    table = search_result.table
    file_name =os.path.join(outdir, "search_resut_tpf_%s.csv" % target)
    table.write(file_name, format ="csv")

def save_tpf(tpfs, search_result,  outdir):
    """ Save tpfs to fits files

        Args:
            tpfs: Arrays for Target Pixel File objects
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            wcs_arr: Arrays of wcs information 

    """
    file_heads = utils.make_output_name_for_tpfs(search_result)
    for (i, tpf) in enumerate(tpfs):
        file_name = os.path.join(outdir, file_heads [i]) + ".fits"
        tpf.to_fits(file_name )


def load_tpf(search_result, outdir):
    """ Load fits files for tpfs

        Args:
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            tpfs: Arrays for Target Pixel File objects

    """    

    file_heads = utils.make_output_name_for_tpfs(search_result)
    tpfs = []
    for (i, file_name) in enumerate(file_heads):
        file_name = os.path.join(outdir, file_heads [i]).replace(" ","")+ ".fits"
        tpf_now = lk.read(file_name)
        tpf_now = tpf_now[tpf_now.quality<15] 
        tpfs.append(tpf_now)
    return tpfs

def load_data_tpf(input_dir, target):
    """ Load tpfs & search_result

        Args:
            input_dir: Directory for output
            target: Target name
        Returns:
            tpfs: Arrays for Target Pixel File objects

    """        
    result_file_name =os.path.join(input_dir, "search_resut_tpf_%s.csv" % target)
    table_result = Table.read(result_file_name, format='csv', fast_reader=False)
    search_result = lk.SearchResult(table = table_result )
    tpfs = load_tpf(search_result, input_dir)
    return tpfs, search_result
    
def compute_centroid(image, x_cen, y_cen, wd):
    """ Compute centroids for images. The region is limited aroun (x_cen, y_cen)

        Args:
            image: 2D Image
            x_cen: Central position in x for computing centroids
            y_cen: Central position in y for computing centroids
            wd: Width of region for computing centroids
        Returns:
            center_x: x centroid
            center_y: y centroid
    """      
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
    center_x = center_x/flux_all
    center_y = center_y/flux_all
    return center_x, center_y
def compute_centroids_for_images(images, wcs_arr, ra_target, dec_target, wd=2):
    """ Compute centroids for images. The region is limited aroun (x_cen, y_cen)

        Args:
            image: 2D Image
            x_cen: Central position in x for computing centroids
            y_cen: Central position in y for computing centroids
            wd: Width of region for computing centroids
        Returns:
            center_x: x centroid
            center_y: y centroid
    """     

    xy_cen_arr = []
    ra_dec_cen_arr = []

    for (j,wcs_now) in enumerate(wcs_arr):
        sky_coord  = SkyCoord(ra_target, dec_target, frame='icrs', unit='deg')
        x_target, y_target =  wcs_arr[j].world_to_pixel(sky_coord)
        print(x_target, y_target)
        cen_x, cen_y = compute_centroid(images[j], x_target, y_target, wd)
        print(cen_x, cen_y)
        ra_cen, dec_cen = wcs_arr[j].wcs_pix2world([cen_x], [cen_y], 0)
        print(ra_cen, dec_cen)
        xy_cen_arr.append([cen_x, cen_y])
        ra_dec_cen_arr.append([ra_cen[0], dec_cen[0]])

    return xy_cen_arr, ra_dec_cen_arr

def make_difference_image_lightkurve(tpfs, period_inputs = None, aperture_mask= "default", filter_length = 7201):
    """ Differential imaging for periodic signals based 

        Args:
            tpfs: Arrays for Target Pixel File objects
            period_inputs: Arrays of periods for differential imaging. If None, we compute periods from lightcurves
            aperture_mask: Aperture_mask
            filter_length (odd int): Length for filter in lk module 
        
        Returns:
            avg_images: Arrys of averaged images
            diff_images: Arrays of differential images
            periods: Periods used for differential imaging
            folded_lcs: Arrays of folded lightcurves using periods 
            folded_bin_lcs: Arrays of binned & folded lightcurves using periods
    """     

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

        min_timestamps, max_timestamps, folded, folded_bin = differential_imaging.make_timestamps_for_differential_imaging(lc_reduced, period)
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

def make_difference_image_for_tpfs(tpfs, period_inputs = None, aperture_mask= "default", filter_length = 7201):
    """ Differential imaging for periodic signals

        Args:
            tpfs: Arrays for Target Pixel File objects
            period_inputs: Arrays of periods for differential imaging. If None, we compute periods from lightcurves
            aperture_mask: Aperture_mask
            filter_length (odd int): Length for filter in lk module 
        
        Returns:
            avg_images: Arrys of averaged images
            diff_images: Arrays of differential images
            periods: Periods used for differential imaging
            folded_lcs: Arrays of folded lightcurves using periods 
            folded_bin_lcs: Arrays of binned & folded lightcurves using periods
    """     

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

        min_timestamps, max_timestamps, folded, folded_bin = differential_imaging.make_timestamps_for_differential_imaging(lc_reduced, period)
        avg_image, diff_image = differential_imaging.make_diff_image(tpf,min_timestamps, max_timestamps, filter_length )

        diff_images.append(diff_image)
        avg_images.append(avg_image)
        folded_lcs.append(folded)
        folded_bin_lcs.append(folded_bin)
        periods.append(period)

    return avg_images, diff_images, periods, folded_lcs, folded_bin_lcs

