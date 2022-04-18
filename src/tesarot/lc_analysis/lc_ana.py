import os 
import lightkurve as lk
from tesarot.utils import utils
from astropy.table import Table

def plot_lc(lc_collection, out_dir ="./"):
    if not os.path.exists(out_dir):
        os.makedirs(out_dif)
    for lc_now in lc_collection:
        lc_now.plot()

def get_periods(lcs, filer_length = 721):
    periods = []
    for lc in lcs:
        lc = reduce_lc_for_periodogram(lc, filer_length = 721)
        pg = lc.to_periodogram(oversample_factor=1)
        period = pg.period_at_max_power
        periods.append(period)
    return periods

def plot_lcs(lc_collection, target,  file_names, out_dir ="./", filer_length = 721):

    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    for (j, lc_now) in enumerate(lc_collection):

        file_name = os.path.join(out_dir, file_names[j] + "_lc_.png")
        lc = reduce_lc_for_periodogram(lc_now, filer_length = 721)
        lc.plot()
        plt.savefig(file_name)
        plt.close()

        file_name = os.path.join(out_dir, file_names[j] + "_pg_.png")
        pg = lc.to_periodogram(oversample_factor=1)
        pg.plot(view='period', scale='log');
        plt.savefig(file_name)
        plt.close()

    for lc_now in lc_collection:
        lc_now.plot()


def get_tess_lcs(target, mission = "TESS", exptime = 120):
    search_result = lk.search_lightcurve(target = target, mission=mission, exptime = exptime)
    lc_collection = search_result.download_all() 
    return lc_collection, search_result

def save_search_result(search_result, target, outdir):
    table = search_result.table
    file_name =os.path.join(outdir, "search_resut_lc_%s.csv" % target)
    table.write(file_name, format ="csv")

def save_lcs(lcs, search_result,  outdir):
    file_heads = utils.make_output_name_for_lcs(search_result)
    for (i, lc) in enumerate(lcs):
        file_name = os.path.join(outdir, file_heads [i]) + ".fits"
        lc.to_fits(file_name )

def load_lcs(search_result, outdir):
    file_heads = utils.make_output_name_for_lcs(search_result)
    lcs = []
    for (i, lc) in enumerate(file_heads):
        file_name = os.path.join(outdir, file_heads [i]).replace(" ","")+ ".fits"
        print(file_name)
        #lk.LightCurve.read(file_name, format='fits',  time_column='date')
        lcs.append(lk.read(file_name))
    return lcs

def load_data(outdir, target):
    result_file_name =os.path.join(outdir, "search_resut_lc_%s.csv" % target)
    table_result = Table.read(result_file_name, format='csv', fast_reader=False)
    search_result = lk.SearchResult(table = table_result )
    print(search_result )
    lcs = load_lcs(search_result, outdir)
    return lcs, search_result

def reduce_lc_for_periodogram(lc, remove_outlier =True, remove_flatten = True, filer_length = 721):
    if remove_outlier:
        lc, mask_out = lc.remove_outliers(return_mask=True)
    if remove_flatten:
        lc = lc.flatten(window_length=filer_length, polyorder=2)
    return lc


