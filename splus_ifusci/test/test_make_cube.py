import getpass

import astropy.units as u

import splusdata
from splus_ifusci import SCube, make_RGB_with_overlay

def make_cube_ngc1087(conn):
    # NGC1087
    galaxy = 'NGC1087'
    coords = ['02:46:25.15', '-00:29:55.45']
    size = 600
    # Connect with S-PLUS
    scube = SCube(galaxy, coords, size, conn=conn,
                  coord_unit=(u.hourangle, u.degree))
    scube.download_stamps(redo=True)
    scube.make_cube(redo=True)
    halpha, halpha_err = scube.calc_halpha()
    # Making RGB image
    flam = scube.get_flam().value
    rgb_bands = ["I", "R", "G"]
    rgb = [flam[scube.bands.index(b)] for b in rgb_bands]
    outimg = f"{galaxy}_RGB.png"
    make_RGB_with_overlay(*rgb, outimg, overlay=halpha.value)

if __name__ == "__main__":
    #Connect with S-PLUS
    username = getpass.getuser() # Change to your S-PLUS username
    password = getpass.getpass(f"Password for {username}:")
    conn = splusdata.connect(username, password)
    make_cube_ngc1087(conn)
