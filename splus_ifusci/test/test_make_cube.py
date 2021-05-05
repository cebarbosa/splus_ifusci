import getpass

import astropy.units as u
import splusdata
from splus_ifusci import SCubeMaker

def make_cube_ngc1087(conn):
    # NGC1087
    galaxy = 'NGC1087'
    coords = ['02:46:25.15', '-00:29:55.45']
    size = 600
    # Connect with S-PLUS
    scube = SCubeMaker(galaxy, coords, size, conn=conn,
                       coord_unit=(u.hourangle, u.degree), redo=True)


if __name__ == "__main__":
    # Connect with S-PLUS
    username = getpass.getuser() # Change to your S-PLUS username
    password = getpass.getpass(f"Password for {username}:")
    conn = splusdata.connect(username, password)
    make_cube_ngc1087(conn)
