import os 
import lightkurve as lk
from tesarot.utils import utils
from astropy.table import Table
import matplotlib.pyplot as plt
from tesarot.lc_analysis import lc_ana


def plot_lc(lc_collection, out_dir ="./"):
    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    for lc_now in lc_collection:
        lc_now.plot()

def plot_lcs(lc_collection, target,  file_names, out_dir ="./", filer_length = 721):

    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    for (j, lc_now) in enumerate(lc_collection):

        file_name = os.path.join(out_dir, file_names[j] + "_lc_.png")
        lc = lc_ana.reduce_lc_for_periodogram(lc_now, filer_length = 721)
        lc.plot()
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()

        file_name = os.path.join(out_dir, file_names[j] + "_pg_.png")
        pg = lc.to_periodogram(oversample_factor=1)
        pg.plot(view='period', scale='log');
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()

    for lc_now in lc_collection:
        lc_now.plot()