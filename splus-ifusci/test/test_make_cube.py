import sys
sys.path.append("..")
import os

import splusdata

from make_cube import SCube

def test1():
    # NGC1087
    galaxy = 'NGC1087'
    coords = ['02:46:25.15', '-00:29:55.45']
    size = 600
    conn = splusdata.connect()
    scube = SCube(galaxy, coords, size)

if __name__ == "__main__":
    test1()
