import numpy as np
import lightkurve as lk
from tesarot.utils import utils
import os
from astropy.table import Table


def get_wcs_from_tpfs(tpfs):
    """ Get wcs information for tpfs.

        Args:
            tpfs: Arrays for Target Pixel File objects
        Returns:
            wcs_arr: Arrays of wcs information 

    """

    wcs_arr = []
    for tpf in tpfs:
        wcs_arr.append(tpf.wcs)
    return wcs_arr