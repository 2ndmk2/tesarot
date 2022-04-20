import os 
import lightkurve as lk
from tesarot.utils import utils
from tesarot.tpf_analysis import wcs_tpf
from tesarot.tpf_analysis import ana_tpf
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np


def plot_tpfs_with_gaia(tpfs, target, out_dir="./", gaia_stars = None,file_names =None, header =""):
    """     
    plot tpfs with gaia stars 
        Args:
            tpfs: Arrays for Target Pixel File objects
            target: Target name
            out_dir: Directory for output
            gaia_stars: dataframe for stars
            file_names: Arrays for files names 
            header: Header for outputfile (normally set to be None)
        Returns:
            None
    """     

    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    wcs_arr = wcs_tpf.get_wcs_from_tpfs(tpfs)

    for (j, tpf) in enumerate(tpfs):
        if file_names is not None:
            file_name = os.path.join(out_dir, "%s_img.png" % (file_names[j]) )
        else:
            file_name = os.path.join(out_dir, "%s_%d_%stpf.png" % (target,j, header) )
        ax = plt.subplot(projection=wcs_arr[j])
        plt.imshow(tpf.flux[0].value)
        if gaia_stars is None:
            plt.savefig(file_name,  bbox_inches="tight")
            plt.close()

        else:
            ra, dec = gaia_stars['ra'], gaia_stars['dec']    
            gaia_mag = gaia_stars['phot_g_mean_mag']
            x, y =  wcs_arr[j].all_world2pix(ra, dec, 0)
            plt.scatter(ra[0], dec[0], transform=ax.get_transform('icrs'),  s = 100, color="b", label="target")
            s = np.maximum((19 - gaia_stars['phot_g_mean_mag'])*10, 0)
            plt.scatter(ra, dec,  s = s, color="r", transform=ax.get_transform('icrs'))
            plt.legend()
            plt.savefig(file_name, bbox_inches="tight" )
            plt.close()

def plot_diffimages_with_gaia(ave_images, diff_images, xy_cen_arr,  target, gaia_stars, periods, file_names, wcs_arr =None, out_dir="./"):
    """     
    plot averaged & differential images with gaia stars 
        Args:
            avg_images: Arrys of averaged images
            diff_images: Arrays of differential images
            cen_x: centroid for x
            cen_y: centroid for y
            target: Target name
            gaia_stars: dataframe for stars
            periods: Arrays of periods for differential imaging. If None, we compute periods from lightcurves
            file_names: Arrays for files names 
            wcs_arr: Arrays of wcs information 
            out_dir: Directory for output
        Returns:
            None
    """ 

    if not os.path.exists(out_dir):
        os.makedirs(out_dif)

    ra, dec = gaia_stars['ra'], gaia_stars['dec']    
    gaia_mag = gaia_stars['phot_g_mean_mag']

    for (j, diff_image) in enumerate(diff_images):
        file_name = os.path.join(out_dir, "%s_diffimg.png" % (file_names[j]) )
        s = np.maximum((19 - gaia_mag)*10, 0)
        x, y =  wcs_arr[j].all_world2pix(ra, dec, 0)
        fig, ax = plt.subplots(1,2, figsize =(12,9))
        ax1 = plt.subplot(121, projection = wcs_arr[j])
        ax2 = plt.subplot(122, projection = wcs_arr[j])
        ax1.imshow(ave_images[j], )
        ax1.scatter(x[0],y[0], s =100, color="b", label="target")
        ax1.scatter(x,y, s =s, color="r")
        ax1.legend()
        ax1.set_title('Average image')
        ax2.imshow(diff_image)
        ax2.set_title('Difference image P:%.4f day' % periods[j].value)
        ax2.scatter(x[0],y[0], s =100, color="b", label="target")
        ax2.scatter(xy_cen_arr[j][0], xy_cen_arr[j][1], s = 40, color="g", label="centroid")
        ax2.scatter(x,y, s =s, color="r", alpha = 0.5)
        ax2.legend()
        plt.savefig(file_name, bbox_inches="tight")
        plt.close()

