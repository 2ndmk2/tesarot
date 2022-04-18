
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut



def data_dl(ra, dec, size):
    cutout_coord = SkyCoord(ra, dec, unit="deg")
    sector_table = Tesscut.get_sectors(coordinates=cutout_coord)
    hdulist = []
    print("sectors:", sector_table["sector"].data)
    for sector_num in sector_table["sector"]:
        print(sector_num)
        hdulist.append(Tesscut.get_cutouts(coordinates=cutout_coord, size=size, sector = sector_num)[0])
    return hdulist


def sector_name_convert(sector_name):
    return int(sector_name.replace("TESS Sector ", ""))