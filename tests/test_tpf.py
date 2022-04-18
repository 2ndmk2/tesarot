from tesarot.lc_analysis import lc_ana
from tesarot.tpf_analysis import ana_tpf
from tesarot.tpf_analysis import plotter_tpf
from tesarot.tpf_analysis import wcs_tpf
from tesarot.utils import utils
import astropy.units as u

import os
import numpy as np

target_for_dl=  "TIC 17198188"
target = target_for_dl.replace(" ","")
out_dir = "./data"
out_dir_for_target = utils.make_output_folder_for_target(target, out_dir)
radius_gaia_search = (12/2.0) * 20*u.arcsec*1.5


## tpf ana
result_file_name_tpf =os.path.join(out_dir_for_target, "search_resut_tpf_%s.csv" % target)
if os.path.exists(result_file_name_tpf):
    tpfs, search_result =ana_tpf.load_data_tpf(out_dir_for_target, target)

else:
    tpfs, search_result = ana_tpf.get_tess_tpf(target_for_dl)
    ana_tpf.save_tpf(tpfs, search_result,  out_dir_for_target)
    ana_tpf.save_search_result_tpf(search_result, target, out_dir_for_target)

ra, dec, gaia_mag, gaia_stars = utils.get_gaia_stars(target, radius_gaia_search)
plotter_tpf.plot_tpfs_with_gaia(tpfs, target, out_dir_for_target, gaia_stars)

## differential imaging 
avg_images, diff_images, periods, folded_lcs, folded_bin_lcs = ana_tpf.make_difference_image_for_tpfs(tpfs, aperture_mask ="default")
wcs_arr = wcs_tpf.get_wcs_from_tpfs(tpfs)
file_heads_tpf = utils.make_output_name_for_tpfs(search_result)
plotter_tpf.plot_diffiamges_with_gaia(avg_images, diff_images, target, \
    gaia_stars, periods,file_heads_tpf,  wcs_arr, out_dir=out_dir_for_target)