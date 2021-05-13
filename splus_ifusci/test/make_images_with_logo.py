import getpass
import os

from PIL import Image
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

import splusdata
from splus_ifusci import SCubeMaker, make_RGB_with_overlay

def make_large_galaxies(conn):
    galaxies = ['NGC4030', "NGC3464", "NGC3314A", "NGC3511", "NGC0428"]
    coordinates = np.array([[180.0986, -1.1002],
                       [163.6666945533915, -21.06532985884086],
                       [159.30363493989583, -27.683956942836808],
                       [165.8485378153153, -23.083694639214354],
                       [18.23199,  0.98148]]) * u.degree
    logo = Image.open("splus_logo.gif")
    for i in range(len(galaxies)):
        galaxy = galaxies[i]
        coords = SkyCoord(*coordinates[i])
        print(coords.ra, coords.dec)
        size = 600
        scube = SCubeMaker(galaxy, coords, size, conn=conn,
                           coord_unit=(u.hourangle, u.degree))
        scube.download_stamps(redo=False)
        scube.make_cube()
        halpha, halpha_err = scube.calc_halpha()
        # Making RGB image
        flam = scube.get_flam().value
        rgb_bands = ["I", "R", "G"]
        rgb = [flam[scube.bands.index(b)] for b in rgb_bands]
        outimg = f"{galaxy}_RGB.jpg"
        make_RGB_with_overlay(*rgb, outimg, overlay=halpha.value)
        img = Image.open(outimg)

        # Including logo
        logo = Image.open("splus_logo.gif").convert("RGBA")
        l, h = logo.size
        ln = int(img.size[0] / 3.)
        hn = int(ln * h / l)
        logon = logo.resize((ln, hn))
        img.paste(logon, (img.size[0]-ln, img.size[1]-hn), logon)
        img.save(outimg.replace(".jpg", "_logo.jpg"))


if __name__ == "__main__":
    #Connect with S-PLUS
    username = getpass.getuser() # Change to your S-PLUS username
    # password = getpass.getpass(f"Password for {username}:")
    # conn = splusdata.connect(username, password)
    conn = None
    make_large_galaxies(conn)
