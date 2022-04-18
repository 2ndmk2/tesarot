import os 
import lightkurve as lk
from tesarot.utils import utils
from astropy.table import Table



def get_tess_lcs(target, mission = "TESS", exptime = 120):
    """ Download lightcurves for a particular target .

        Args:
            target: Target name
            mission: Mission name
            exptime: Cadence for data
        Returns:
            lc_collection: Arrays of Lightcurve objects 
            search_result: Table containing search result for data 

    """

    search_result = lk.search_lightcurve(target = target, mission=mission, exptime = exptime)
    lc_collection = search_result.download_all() 
    return lc_collection, search_result

def save_search_result(search_result, target, outdir):
    """ Save "search result" named as "file_name".

        Args:
            search_result: Table containing search result for data 
            target: Target name
            out_dir: Directory for output
        Returns:
            None

    """    

    table = search_result.table
    file_name =os.path.join(outdir, "search_resut_lc_%s.csv" % target)
    table.write(file_name, format ="csv")

def save_lcs(lcs, search_result,  outdir):
    """ Save Lightcurve objects to fits files.

        Args:
            lcs: Arrays of Lightcurve object
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            None

    """  

    file_heads = utils.make_output_name_for_lcs(search_result)
    for (i, lc) in enumerate(lcs):
        file_name = os.path.join(outdir, file_heads [i]) + ".fits"
        lc.to_fits(file_name )

def load_lcs(search_result, outdir):
    """ Load fits files for Lightcurve objects.

        Args:
            search_result: Table containing search result for data 
            out_dir: Directory for output
        Returns:
            lcs: Arrays of Lightcurve object

    """  
    file_heads = utils.make_output_name_for_lcs(search_result)
    lcs = []
    for (i, lc) in enumerate(file_heads):
        file_name = os.path.join(outdir, file_heads [i]).replace(" ","")+ ".fits"
        print(file_name)
        #lk.LightCurve.read(file_name, format='fits',  time_column='date')
        lcs.append(lk.read(file_name))
    return lcs

def load_data(outdir, target):
    """ Load fits files for Lightcurve objects.

        Args:
            out_dir: Directory for output
            target: Target name
        Returns:
            lcs: Arrays of Lightcurve object
            search_result: Table containing search result for data 
    """  

    result_file_name =os.path.join(outdir, "search_resut_lc_%s.csv" % target)
    table_result = Table.read(result_file_name, format='csv', fast_reader=False)
    search_result = lk.SearchResult(table = table_result )
    lcs = load_lcs(search_result, outdir)
    return lcs, search_result


def get_periods(lcs, filter_length = 721):
    """ Get arrays for periods for Lightcurve objects.

        Args:
            lcs: Arrays of Lightcurve object
            filter_length (int): Length for filter in lk module 
        Returns:
            periods: Arrays of periods (u.day)
    """  

    periods = []
    for lc in lcs:
        lc = reduce_lc_for_periodogram(lc, filter_length = 721)
        pg = lc.to_periodogram(oversample_factor=1)
        period = pg.period_at_max_power
        periods.append(period)
    return periods

def reduce_lc_for_periodogram(lc, remove_outlier =True, remove_flatten = True, filter_length = 721):
    """ Remove outliers & flatten lcs. 

        Args:
            lc: Lightcurve object
            remove_outlier: If True, we remove outliers
            remove_flatten: If True, we flatten lightcurve
            filter_length (int): Length for filter in lk module 
        Returns:
            periods: Arrays of periods (u.day)
    """  

    if remove_outlier:
        lc, mask_out = lc.remove_outliers(return_mask=True)
    if remove_flatten:
        lc = lc.flatten(window_length=filter_length, polyorder=2)
    return lc


