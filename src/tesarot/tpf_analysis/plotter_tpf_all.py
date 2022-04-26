import os 
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from tesarot.lc_analysis import lc_ana
import matplotlib.gridspec as gridspec


def make_output_name_for_ffis_periodogram(sector_info, target_name, out_dir):

    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            files: Arrays for files names of tpfs
    """ 

    sectors = sector_info["sectorName"]
    file_names = []

    for i in range(len(sectors)):
        file_name = "ffi_TESS_%s_%s_period_all.pdf" %(target_name, sectors[i])
        file_names.append(os.path.join(out_dir, file_name))

    return file_names


def make_output_name_for_tpfs_periodogram(search_result, out_dir):

    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            files: Arrays for files names of tpfs
    """ 

    missions = search_result.mission
    targets = search_result.target_name
    exp = search_result.exptime
    files = []
    for i in range(len(missions)):
        mission_now = missions[i].replace(" ", "_")
        exp_now = str(exp[i]).replace("s", "")
        file_name = "tpf_%s_%s_%ds_period_all.pdf" %(mission_now , targets[i], int(float(exp_now)))
        files.append(os.path.join(out_dir, file_name))
    return files


def make_output_name_for_tpfs_lc(search_result, out_dir):

    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            files: Arrays for files names of tpfs
    """ 

    missions = search_result.mission
    targets = search_result.target_name
    exp = search_result.exptime
    files = []
    for i in range(len(missions)):
        mission_now = missions[i].replace(" ", "_")
        exp_now = str(exp[i]).replace("s", "")
        file_name = "tpf_%s_%s_%ds_lc_all.pdf" %(mission_now , targets[i], int(float(exp_now)))
        files.append(os.path.join(out_dir, file_name))
    return files


def periodogram_for_tpf(tpf, period_ref = None, window_length = 6001, period_min = 0.1, \
    period_max = 10, file_name ="test.pdf", aperture_now = None):
    """ Plotter for peridogram for tpf

    Params:
        tpf: Target Pixel File object
        period_ref: Reference period in plot
        window_length (odd int): Length for filter in lk module 
        period_min: Minimum period in plot
        period_max: Max period in plot
        file_name: output file_name
        aperture_now: aperture_region for target
    Returns:
        None

    """

    nx, ny = np.shape(tpf.flux[0])    
    gs = gridspec.GridSpec(
        nx, ny, wspace=0.01, hspace=0.01
    )

    fig = plt.figure()
    ax = plt.gca()
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if aperture_now is None:
        aperture_now = tpf.pipeline_mask

    for k in range(nx * ny):        
        x, y = np.unravel_index(k, (nx, ny))
        aperture_mask_now = np.zeros((nx, ny))
        aperture_mask_now[x,y] = 1
        lc = tpf.to_lightcurve(aperture_mask=aperture_mask_now )
        pg_raw = lc.to_periodogram(oversample_factor=1)
        lc = lc_ana.reduce_lc_for_periodogram(lc)
        lc = lc.remove_nans().flatten(window_length = window_length).remove_outliers()  
        pg = lc.to_periodogram(oversample_factor=1)

        if aperture_now[x, y]:
            rc = {"axes.linewidth": 2, "axes.edgecolor": "red"}
        else:
            rc = {"axes.linewidth": 1}
            

        with plt.rc_context(rc=rc):
            gax = fig.add_subplot(gs[nx - x - 1, y])

        gax.plot(
            pg.period.value,
            pg.power.value/np.nanmax(pg.power.value),
            marker="None",
            color="k"
        )    
        gax.plot(
            pg_raw.period.value,
            pg_raw.power.value/np.nanmax(pg_raw.power.value),
            marker="None",
            color="r",
            alpha = 0.5
        )    

        if period_ref is not None:
            gax.axvline(x=period_ref, color="r", lw = 1, linestyle="dashed", zorder = -100)
        gax.set_xscale("log")
        gax.set_xlim( period_min ,  period_max)
        gax.set_xticklabels("")
        gax.set_yticklabels("")
        gax.set_xticks([])
        gax.set_yticks([])    
    fig.set_size_inches((15, 15))
    plt.savefig(file_name,  bbox_inches="tight")

def periodogram_for_tpfs(tpfs, file_names, period_refs = None, window_length = 6001, period_min = 0.1, \
        period_max = 10, aperture_now = None):
    """ Plotter for peridogram for tpfs

    Params:
        tpfs: Array of Target Pixel File object
        file_names: Array of output file_names
        window_length (odd int): Length for filter in lk moadule 
        aperture_now: aperture_region for target

    Returns:
        None

    """

    for (i, tpf) in enumerate(tpfs):

        if period_refs is None:
            periodogram_for_tpf(tpf, period_ref = None, window_length =window_length, 
                period_min = period_min, period_max  = period_max , file_name = file_names[i], aperture_now = aperture_now)
        else:
            periodogram_for_tpf(tpf, period_ref = period_refs[i].value, window_length =window_length, 
                period_min = period_min, period_max  = period_max , file_name = file_names[i], aperture_now = aperture_now) 


def lcs_for_tpf(tpf, window_length = 6001, file_name ="test.pdf", aperture_now=None):
    """ Plotter for lightcurves for tpf

    Params:
        tpf: Target Pixel File object
        window_length (odd int): Length for filter in lk module 
        file_name: output file_name
    Returns:
        None

    """

    nx, ny = np.shape(tpf.flux[0])    
    gs = gridspec.GridSpec(
        nx, ny, wspace=0.01, hspace=0.01
    )

    fig = plt.figure()
    ax = plt.gca()
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    if aperture_now is None:
        aperture_now = tpf.pipeline_mask

    for k in range(nx * ny):        
        x, y = np.unravel_index(k, (nx, ny))
        aperture_mask_now = np.zeros((nx, ny))
        aperture_mask_now[x,y] = 1
        lc_raw = tpf.to_lightcurve(aperture_mask=aperture_mask_now )
        lc = lc_ana.reduce_lc_for_periodogram(lc_raw)
        lc = lc.remove_nans().flatten(window_length = window_length).remove_outliers()  

        if aperture_now[x, y]:
            rc = {"axes.linewidth": 2, "axes.edgecolor": "red"}
        else:
            rc = {"axes.linewidth": 1}

        with plt.rc_context(rc=rc):
            gax = fig.add_subplot(gs[nx - x - 1, y])

        gax.plot(
            lc_raw.time.value,
            lc_raw.flux.value/np.nanmax(lc_raw.flux.value),
            marker="None",
            color="k"
        )    
        gax.plot(
            lc.time.value,
            lc.flux.value/np.nanmax(lc.flux.value),
            marker="None",
            color="r",
            alpha = 0.5
        )          

        gax.set_xticklabels("")
        gax.set_yticklabels("")
        gax.set_xticks([])
        gax.set_yticks([])    
    fig.set_size_inches((15, 15))
    plt.savefig(file_name,  bbox_inches="tight")


def lcs_for_tpfs(tpfs, file_names,window_length = 6001, aperture_now=None):
    """ Plotter for lightcurves for tpfs

    Params:
        tpfs: Array of Target Pixel File object
        file_names: Array of output file_names
        window_length (odd int): Length for filter in lk moadule 
        aperture_now: aperture_region for target
    Returns:
        None

    """

    for (i, tpf) in enumerate(tpfs):

            lcs_for_tpf(tpf, window_length , file_names[i], aperture_now) 



