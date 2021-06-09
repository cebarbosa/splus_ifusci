import getpass
import os

from PIL import Image
from tqdm import tqdm
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

import splusdata
from splus_ifusci import SCube, make_RGB_with_overlay

def make_large_galaxies(conn):
    galaxies = ["NGC1365", "FCC37", "ARP244", "quasar_3C273", 'NGC4030',
                "NGC3464",
                "NGC3314A",
               "NGC3511",
                "NGC0428", "NGC7089", "HydraCluster"]
    coordinates = np.array([[53.40166, -36.140277],
                            [51.33458, -36.385],
                            [180.47208, -18.87694],
                            [187.2779, 2.0525],
                        [180.0986, -1.1002],
                       [163.6666945533915, -21.06532985884086],
                       [159.30363493989583, -27.683956942836808],
                       [165.8485378153153, -23.083694639214354],
                       [18.23199,  0.98148],
                       [323.3625, -0.823333],
                       [159.17416, -27.525444]]) * u.degree
    sizes = [600] * len(galaxies)
    sizes[0] = 1800
    sizes[1] = 900
    sizes[2] = 2400
    sizes[-2] = 2400
    sizes[-1] = 2400
    for i in tqdm(range(len(galaxies)), desc="Processing objects"):
        galaxy = galaxies[i]
        coords = SkyCoord(*coordinates[i])
        size = sizes[i]
        scube = SCube(galaxy, coords, size, conn=conn,
                      coord_unit=(u.hourangle, u.degree))
        scube.download_stamps(redo=False)
        scube.make_cube()
        halpha, halpha_err = scube.calc_halpha()
        # Making RGB image
        flam = scube.get_flam().value
        rgb_bands = ["I", "R", "G"]
        rgb = [flam[scube.bands.index(b)] for b in rgb_bands]
        outimg = f"{galaxy}_{size}x{size}_RGB.jpg"
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
    password = getpass.getpass(f"Password for {username}:")
    conn = splusdata.connect(username, password)
    # conn = None
    make_large_galaxies(conn)
