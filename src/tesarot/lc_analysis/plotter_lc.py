import os 
import lightkurve as lk
from tesarot.utils import utils
from astropy.table import Table
import matplotlib.pyplot as plt
from tesarot.lc_analysis import lc_ana


def plot_lc(lc_collection):
    """ Plot lightcurves.

        Args:
            lc_colletion: Arrays of Lightcurve object
            out_dir: Directory for output
        Returns:
            None

    """


    for lc_now in lc_collection:
        lc_now.plot()

def plot_lcs(lc_collection, target,  file_names, out_dir ="./", filter_length = 7201, exptime = 120):

    """ Plot lightcurves and their periodgrams.

        Args:
            lc_colletion: Arrays of Lightcurve object
            target: Name of target
            file_names: Arrays of names for outputfiles
            out_dir: Output directory
            filter_length (int): Length for filter in lk module 

        Returns:
            None

    """
    filter_width_day = filter_length * exptime/(3600*24.0)

    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    for (j, lc_now) in enumerate(lc_collection):

        file_name = os.path.join(out_dir, file_names[j] + "_lc_raw.png")
        title_now = "raw lc %s" % target
        lc_now.plot(title = title_now)
        plt.savefig(file_name, bbox_inches="tight", title = title_now)
        plt.close()

        file_name = os.path.join(out_dir, file_names[j] + "_pg_raw.png")
        pg = lc_now.to_periodogram(oversample_factor=1)
        title_now = "reduced periodogram %s" % target
        pg.plot(view='period', scale='log', title = title_now);
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()        

        file_name = os.path.join(out_dir, file_names[j] + "_lc_.png")
        lc = lc_ana.reduce_lc_for_periodogram(lc_now, filter_length = filter_length)
        title_now = "reduced lc %s:(High pass %f .days)" % (target, filter_width_day)
        lc.plot( title = title_now)
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()

        file_name = os.path.join(out_dir, file_names[j] + "_pg_.png")
        pg = lc.to_periodogram(oversample_factor=1)
        title_now = "reduced periodogram %s:(High pass %fdays)" % (target, filter_width_day)
        pg.plot(view='period', scale='log', title = title_now );
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()

    for lc_now in lc_collection:
        lc_now.plot()
