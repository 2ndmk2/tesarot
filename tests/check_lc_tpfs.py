from tesarot.lc_analysis import lc_ana
from tesarot.tpf_analysis import ana_tpf
from tesarot.lc_analysis import plotter_lc
from tesarot.tpf_analysis import plotter_tpf
from tesarot.tpf_analysis import plotter_tpf_all
from tesarot.tpf_analysis import wcs_tpf
from tesarot.utils import utils
import astropy.units as u
import os


def main(target_for_dl, exptime=120, out_dir ="./data",  filter_length =7201):

    target = target_for_dl.replace(" ","")
    out_dir_for_target = utils.make_output_folder_for_target(target, out_dir)

    ## lc analysis
    print("lc analysis")
    result_file_name =os.path.join(out_dir_for_target, "search_result_lc_%s.csv" % target)
    if os.path.exists(result_file_name):
        lc_collection, search_result = lc_ana.load_data(out_dir_for_target, target)
    else:
        lc_collection, search_result = lc_ana.get_tess_lcs(target_for_dl, exptime=exptime)
        if lc_collection is None:
            print("no data:",target_for_dl)
            return None
        lc_ana.save_lcs(lc_collection, search_result, out_dir_for_target)
        lc_ana.save_search_result(search_result, target, out_dir_for_target)
    periods = lc_ana.get_periods(lc_collection, filter_length )


    ## tpf analysis
    print("tpf analysis")

    result_file_name_tpf =os.path.join(out_dir_for_target, "search_result_tpf_%s.csv" % target)
    if os.path.exists(result_file_name_tpf):
        tpfs, search_result_tpf =ana_tpf.load_data_tpf(out_dir_for_target, target)

    else:
        tpfs, search_result_tpf = ana_tpf.get_tess_tpf(target_for_dl, exptime=exptime)
        if tpfs is None:
            print("no data:",target_for_dl)
            return None        
        ana_tpf.save_tpf(tpfs, search_result_tpf,  out_dir_for_target)
        ana_tpf.save_search_result_tpf(search_result_tpf, target, out_dir_for_target)

    files = plotter_tpf_all.make_output_name_for_tpfs_lc(search_result_tpf, out_dir_for_target )
    plotter_tpf_all.lcs_for_tpfs(tpfs, files,  filter_length )

    files = plotter_tpf_all.make_output_name_for_tpfs_periodogram(search_result_tpf, out_dir_for_target)
    plotter_tpf_all.periodogram_for_tpfs(tpfs, files, periods, filter_length)

if __name__ == "__main__":

    targets_for_dl=  ["TIC 121462259", "TIC 130414928", "TIC 141522677", "TIC 144069454"]
    exptime=120 # seconds
    out_dir ="./data"
    filter_length_day = 5.0 #day
    filter_length = utils.convert_filter_width(filter_length_day, exptime)
    print("filter_wd:", filter_length)

    if filter_length %2==0:
        print("filter_length must be odd number!!")
    else:
        main(target_for_dl, exptime, out_dir, filter_length )
