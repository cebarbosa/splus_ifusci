import getpass

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

import splusdata
from splus_ifusci import SCubeMaker, SPLUSEmLine, make_RGB_with_overlay

def make_cube_ngc1087(conn):
    # NGC1087
    galaxy = 'NGC1087'
    coords = ['02:46:25.15', '-00:29:55.45']
    size = 600
    # Connect with S-PLUS
    scube = SCubeMaker(galaxy, coords, size, conn=conn,
                       coord_unit=(u.hourangle, u.degree), redo=False)
    flam = scube.get_flam().value
    # Calculating H-alpha
    wline = 6562.8 * u.AA
    bands = ["R", "F660", "I"]
    halpha_estimator = SPLUSEmLine(wline, bands)
    idx = np.array([scube.bands.index(band) for band in scube.bands])
    flux_halpha = halpha_estimator.flux_3F(flam[idx])
    vmin = np.nanpercentile(flux_halpha.value, 10)
    vmax = np.percentile(flux_halpha.value, 99)
    # Making RGB image
    rgb_bands = ["I", "R", "G"]
    rgb = [flam[scube.bands.index(b)] for b in rgb_bands]
    outimg = f"{galaxy}_RGB.png"
    make_RGB_with_overlay(*rgb, outimg)
    # plt.imshow(flux_halpha.value, origin="lower", vmin=vmin, vmax=vmax)
    # plt.colorbar()
    # plt.show()




if __name__ == "__main__":
    # Connect with S-PLUS
    # username = getpass.getuser() # Change to your S-PLUS username
    # password = getpass.getpass(f"Password for {username}:")
    # conn = splusdata.connect(username, password)
    conn = None
    make_cube_ngc1087(conn)
