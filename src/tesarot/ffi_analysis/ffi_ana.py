import os 
import lightkurve as lk
from tesarot.utils import utils
from astropy.table import Table
from astroquery.mast import Catalogs
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut


def get_ra_dec(target):
    """ Obtain ra & dec for object

        Args:
            target: Target name
        Returns:
            ra: right ascention
            dec: declination
    """  

    catalogTIC = Catalogs.query_object("TIC 121638493", radius=0.01 * u.arcsec, catalog="TIC")
    ra = catalogTIC["ra"]
    dec = catalogTIC["dec"]
    return ra, dec

def load_data(outdir, target):
    """ Load fits files for Pixel objects. 

        Args:
            out_dir: Directory for output
            target: Target name
        Returns:
            tpfs: Arrays of Pixel objects
            sector_info: Table containing search result for data 
    """  

    result_file_name =os.path.join(outdir, "search_resut_ffi_%s.csv" % target)
    sector_info = Table.read(result_file_name, format='csv', fast_reader=False)
    file_names = make_output_name_for_ffi(sector_info, target)
    tpfs = ffi_load(file_names, outdir)

    return tpfs, sector_info


def load_sector_info(outdir, target):
    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_info: Table containing sector information

        Returns:
            None

    """
    file_name =os.path.join(outdir, "search_resut_ffi_%s.csv" % target)
    sector_table.write(file_name, format ="csv")

def save_sector_info(sector_info, out_dir, target):
    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_info: Table containing sector information
        Returns:
            None

    """
    file_name =os.path.join(out_dir, "search_resut_ffi_%s.csv" % target)
    sector_info.write(file_name, format ="csv")

def ffi_dl(ra, dec, size):
    """ Download FFI 
        
        Args:
            ra: right ascention
            dec: declination

        Returns:
            size: size of FFI (side)

    """


    cutout_coord = SkyCoord(ra, dec, unit="deg")
    sector_info = Tesscut.get_sectors(coordinates=cutout_coord)
    hdulist = []
    for sector_num in sector_info["sector"]:
        hdulist.append(Tesscut.get_cutouts(coordinates=cutout_coord, size=size, sector = sector_num)[0])
    return hdulist, sector_info

def make_output_name_for_ffi(sector_info, target_name):

    """ Make file names for target-pixel files (tpfs) based on search result
        
        Args:
            search_info: Table containing sector information
            target_name: Target name
        Returns:
            file_names: Arrays for files names of ffis
    """ 

    sectors = sector_info["sectorName"]
    file_names = []

    for i in range(len(sectors)):
        file_name = "ffi_TESS_%s_%s" %(target_name, sectors[i])
        file_names.append(file_name)

    return file_names

def ffi_save(hdulist, file_names, out_dir):
    """ Save FFI files
        
        Args:
            hdul_list: Arrays for FFI data
            file_names: Arrays for files names of fits
            out_dir: Directory for output


    """ 

    for i in range(len(file_names)):
        tpf_name = os.path.join(out_dir, "%s.fits" % (file_names[i]))

        if not os.path.exists(tpf_name):
            hdulist[i].writeto(tpf_name)

    return None

def ffi_load(file_names, out_dir):
    """ Load FFI files
        
        Args:
            file_names: Arrays for files names of fits
            out_dir: Directory for output

        Returns:
            tpfs: Arrays of Pixel objects
    """ 
    tpfs = []

    for i in range(len(file_names)):

        tpf_name = os.path.join(out_dir, "%s.fits" % (file_names[i]))
        tpfs.append( lk.TessTargetPixelFile(tpf_name))
    return tpfs 
