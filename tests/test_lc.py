from tesarot.lc_analysis import lc_ana
from tesarot.tpf_analysis import ana_tpf
from tesarot.lc_analysis import plotter_lc
from tesarot.utils import utils
import astropy.units as u

import os

target_for_dl=  "TIC 17198188"
target = target_for_dl.replace(" ","")
out_dir = "./data"
out_dir_for_target = utils.make_output_folder_for_target(target, out_dir)
radius_gaia_search = (12/2.0) * 20*u.arcsec*1.5

## lc ana
result_file_name =os.path.join(out_dir_for_target, "search_resut_lc_%s.csv" % target)
if os.path.exists(result_file_name):
    lc_collection, search_result = lc_ana.load_data(out_dir_for_target, target)
else:
    lc_collection, search_result = lc_ana.get_tess_lcs(target_for_dl)
    lc_ana.save_lcs(lc_collection, search_result, out_dir_for_target)
    lc_ana.save_search_result(search_result, target, out_dir_for_target)
plotter_lc.plot_lcs(lc_collection, target, out_dir_for_target)
