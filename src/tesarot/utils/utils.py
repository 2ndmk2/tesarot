
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
import os
from astroquery.mast import Catalogs
import astropy.units as u

def get_gaia_stars(target, radius):
    stars = Catalogs.query_object(target, radius=radius, catalog='Gaia')    
    ra, dec = stars['ra'], stars['dec']    
    gaia_mag = stars['phot_g_mean_mag']
    return ra, dec, gaia_mag, stars


def data_dl(ra, dec, size):
    cutout_coord = SkyCoord(ra, dec, unit="deg")
    sector_table = Tesscut.get_sectors(coordinates=cutout_coord)
    hdulist = []
    print("sectors:", sector_table["sector"].data)
    for sector_num in sector_table["sector"]:
        print(sector_num)
        hdulist.append(Tesscut.get_cutouts(coordinates=cutout_coord, size=size, sector = sector_num)[0])
    return hdulist

def make_output_folder_for_target(target, outdir_root):

    outdir_for_tareget = os.path.join(outdir_root, target)
    if not os.path.exists(outdir_for_tareget):
        os.makedirs(outdir_for_tareget)
    return outdir_for_tareget

def make_output_name_for_lcs(search_result):
    missions = search_result.mission
    targets = search_result.target_name
    exp = search_result.exptime
    files = []
    for i in range(len(missions)):
        mission_now = missions[i].replace(" ", "_")
        exp_now = str(exp[i]).replace(" s", "s")
        file_name = "lc_%s_%s_%s" %(mission_now , targets[i], exp_now )
        files.append(file_name)
    return files

def make_output_name_for_tpfs(search_result):
    missions = search_result.mission
    targets = search_result.target_name
    exp = search_result.exptime
    files = []
    for i in range(len(missions)):
        mission_now = missions[i].replace(" ", "_")
        exp_now = str(exp[i]).replace(" s", "s")
        file_name = "tpf_%s_%s_%s" %(mission_now , targets[i], exp_now )
        files.append(file_name)
    return files

def sector_name_convert(sector_name):
    return int(sector_name.replace("TESS Sector ", ""))