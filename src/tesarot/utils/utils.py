
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
import os
from astroquery.mast import Catalogs
import astropy.units as u
import numpy as np

def get_gaia_stars(target, radius):

    """ Get nearby stars from Gaia 
        Args:
            target: Target name
            radius: Search radius (give it with unit of "u.arcsec" in Astropy)
        Returns:
            ra: Arrays of ra
            dec: Arrays of dec
            gaia_mag: Arrays of Gaia magnitude
            stars: dataframe for stars
    """  

    stars = Catalogs.query_object(target, radius=radius, catalog='Gaia')    
    ra, dec = stars['ra'], stars['dec']    
    gaia_mag = stars['phot_g_mean_mag']
    return ra, dec, gaia_mag, stars

def data_dl(ra, dec, size):

    """ Get FFI data
        Args:
            ra: ra
            dec: dec
            size: cutout size (pixel)
        Returns:
            hdulist: arrays for FFI frames
    """  

    cutout_coord = SkyCoord(ra, dec, unit="deg")
    sector_table = Tesscut.get_sectors(coordinates=cutout_coord)
    hdulist = []
    print("sectors:", sector_table["sector"].data)
    for sector_num in sector_table["sector"]:
        print(sector_num)
        hdulist.append(Tesscut.get_cutouts(coordinates=cutout_coord, size=size, sector = sector_num)[0])
    return hdulist

def make_output_folder_for_target(target, outdir_root):

    """ Make directories for target folders
        Args:
            target: Target name
            outdir_root: Path to output directory 
        Returns:
            outdir_for_target: Path to output directory
    """      

    outdir_for_target = os.path.join(outdir_root, target)
    if not os.path.exists(outdir_for_target):
        os.makedirs(outdir_for_target)
    return outdir_for_target

def make_output_name_for_lcs(search_result):

    """ Make file names for lightcurves based on search result
        Args:
            search_result: Table containing search result for data 
        Returns:
            files: Arrays for files names of lcs
    """  

    missions = search_result.mission
    targets = search_result.target_name
    exp = search_result.exptime
    files = []
    for i in range(len(missions)):
        mission_now = missions[i].replace(" ", "_")
        exp_now = str(exp[i]).replace("s", "")
        file_name = "lc_%s_%s_%ds" %(mission_now , targets[i], int(float(exp_now)) )
        files.append(file_name)
    return files

def make_output_name_for_tpfs(search_result):

    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_result: Table containing search result for data 
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
        file_name = "tpf_%s_%s_%ds" %(mission_now , targets[i], int(float(exp_now)))
        files.append(file_name)
    return files

def sector_name_convert(sector_name):
    """ Remove TESS Sector from str

        Args:
            sector_name: str for sector name
        Returns:
            sector_name_changed: changed name
    """
    sector_name_changed = int(sector_name.replace("TESS Sector ", ""))
    return sector_name_changed

def save_results(gaia_stars, periods, ra_dec_cen_arr, out_dir):

    file_name = os.path.join(out_dir, "summary")
    np.savez(file_name, gaia_ra=gaia_stars["ra"],  gaia_dec=gaia_stars["dec"], gaig_mag = gaia_stars["phot_g_mean_mag"], \
            period_arr = [period.value for period in periods], ra_dec_cen_arr = ra_dec_cen_arr)
def convert_filter_width(filter_length_day, exptime):
    """ Get filter window for FFT in unit of cadence for that in unit of day:
        
        Args:
            filter_day: Filter window for FFT in unit of day
            exptime: cadence in unit of second:
        Retruns:
            filter_width: Filter window for FFT in unit of cadence
    """

    filter_length = int(filter_length_day * 3600*24/exptime)
    if filter_length %2 ==0:
        filter_length  += 1
    return filter_length

def load_results(out_dir):

    file_name = os.path.join(out_dir, "summary.npz")
    result = np.load(file_name)
    return result