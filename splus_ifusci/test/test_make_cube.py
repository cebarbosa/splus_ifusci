import sys
sys.path.append("..")
import os

import astropy.units as u
import splusdata
from splus_ifusci import SCubeMaker

from credentials import user, password

def test1():
    galaxy = 'NGC1374'
    coords = ['03:35:16.598', '-35:13:34.50'] # real coordinates
    size = 900
    # Connect with S-PLUS
    conn = splusdata.connect(user, password)
    scube = SCubeMaker(galaxy, coords, size, conn=conn)

def test2():
    # NGC1087
    galaxy = 'NGC1087'
    coords = ['02:46:25.15', '-00:29:55.45']
    size = 600
    # Connect with S-PLUS
    conn = splusdata.connect(user, password)
    scube = SCubeMaker(galaxy, coords, size, conn=conn,
                       coord_unit=(u.hourangle, u.degree), redo=True)


if __name__ == "__main__":
    test2()
