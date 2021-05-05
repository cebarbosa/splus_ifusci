import os
import warnings

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from scipy.interpolate import RectBivariateSpline

warnings.simplefilter('ignore', category=AstropyWarning)
np.seterr(divide='ignore', invalid='ignore')

class SCubeMaker():
    def __init__(self, obj, coords, size, conn=None, wdir=None, zpref=None,
                 redo=False, bands=None, verbose=True, coord_unit=None):
        """ Produces S-PLUS data cubes directly from the project's database.

        Prameters
        ---------
        obj: str
            Identification of the object
        coords: astropy.coordinates.SkyCoord or list
            Coordinates of the object. If list is given, it assumes units in
            degrees for both coordinates.
        size: astropy.Quantity or flot
            Size of the cube.
        conn: splusdata.conn
            Connection with S-PLUS data base using splusdata
        wdir: str
            Working directory
        zpref: str
            Calibration reference. Options are idr3_n4 (latest calibration,
            default) and idr3
        redo: bool
            Produces cube again even if it already existis in the disk.
            Default is False.
        bands: list
            Names of the filters to be included in the cube. Options include
            'U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
            'F861', and 'Z'. Defalt is all bands.
        coord_unit: astropy.unit
            Units to be used in coordinates. Default is degree
        """
        self.obj = obj
        self.coord_unit = u.degree if coord_unit is None else coord_unit
        if isinstance(coords, SkyCoord):
            self.coords = coords
        elif isinstance(coords, list):
                self.coords = SkyCoord(coords[0], coords[1],
                                       unit=self.coord_unit)
        else:
            raise ValueError("Input coordinates should be list or SkyCoord")
        self.size = size
        if isinstance(self.size, u.Quantity):
            self.size_unit = self.size.unit
        else:
            self.size_unit = u.pix
        self.conn = conn
        # General definitions
        self.ps = 0.55 * u.arcsec / u.pixel
        if isinstance(self.size, u.Quantity):
            self.cutsize = int((self.size / self.ps).value)
        else:
            self.cutsize = int(self.size)
        all_bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R',
                      'F660', 'I', 'F861', 'Z']
        if bands is not None:
            assert all([x in all_bands for x in bands]), \
                   "Input bands list should include only S-PLUS filters"
        self.bands = all_bands if bands is None else bands
        self.bands_names = {'U': "$u$", 'F378': "$J378$", 'F395': "$J395$",
                            'F410': "$J410$", 'F430': "$J430$", 'G': "$g$",
                            'F515': "$J515$", 'R': "$r$", 'F660': "$J660$",
                            'I': "$i$", 'F861': "$J861$", 'Z': "$z$"}
        self.wave_eff = {"F378": 3770.0, "F395": 3940.0, "F410": 4094.0,
                         "F430": 4292.0, "F515": 5133.0, "F660": 6614.0,
                         "F861": 8611.0, "G": 4751.0, "I": 7690.0, "R": 6258.0,
                         "U": 3536.0, "Z": 8831.0}
        self.wave = np.array([self.wave_eff[band] for band in
                              self.bands]) *u.Angstrom
        self.flam_unit = u.erg / u.cm / u.cm / u.s / u.AA
        self.fnu_unit = u.erg / u.s / u.cm / u.cm / u.Hz
        # Setting zero point calibration configurations
        self._path = os.path.dirname(os.path.abspath(__file__))
        self.zpref = "idr3_n4" if zpref is None else zpref
        self.zps = self.get_zps()
        self.zpcorr = self.get_zp_correction()
        # Setting directory for stamps
        self.wdir = os.getcwd() if wdir is None else wdir
        self.cutouts_dir = os.path.join(self.wdir, "cutouts")
        if not os.path.exists(self.cutouts_dir):
            os.mkdir(self.cutouts_dir)
        # Producing stamps and cube
        s = int(self.size)
        self.cutnames = [f"{self.obj}_{band}_{s}x{s}{self.size_unit}.fz" for
                         band in self.bands]
        self.wcutnames= [cut.replace(".fz", "_weight.fz") for cut in
                         self.cutnames]
        self.cubename = os.path.join(self.wdir, f"{self.obj}_{s}x{s}"
                                                f"{self.size_unit}.fits")
        status = self.check_infoot()
        if os.path.exists(self.cubename) and not redo:
            return
        if not status:
            if verbose:
                print(f"Data not available for {self.obj}, cube will not be " \
                      f"processed.")
            return
        self.download_stamps()
        self.make_cube()

    def check_infoot(self):
        """ Check if proposed objects are found in the S-PLUS footprint. """
        ra = self.coords.ra.to(u.degree).value
        dec = self.coords.dec.to(u.degree).value
        try:
            hdu = self.conn.get_cut(ra, dec, self.cutsize, self.bands[0])
            return 1
        except:
            return 0

    def download_stamps(self):
        if self.conn is None:
            raise ValueError("A connection with S-PLUS database should be "
                             "provided for missing stamps.")
        for imglist in [self.cutnames, self.wcutnames]:
            for band, stamp in zip(self.bands, imglist):
                outfile = os.path.join(self.cutouts_dir, stamp)
                if os.path.exists(outfile):
                    continue
                ra = self.coords.ra.to(u.degree).value
                dec = self.coords.dec.to(u.degree).value
                hdu = self.conn.get_cut(ra, dec, self.cutsize, band)
                hdu.writeto(outfile)

    def get_zps(self):
        """ Load all tables with zero points for iDR3. """

        zp_dir = os.path.join(self._path, f"assets/{self.zpref}/zps")
        tables = []
        for fname in os.listdir(zp_dir):
            filename = os.path.join(zp_dir, fname)
            data = np.genfromtxt(filename, dtype=None, encoding=None)
            with open(filename) as f:
                h = f.readline().replace("#", "").replace("SPLUS_", "").split()
            table = Table(data, names=h)
            tables.append(table)
        zptable = vstack(tables)
        return zptable

    def get_zp_correction(self):
        """ Get corrections of zero points for location in the field. """
        x0, x1, nbins = 0, 9200, 16
        xgrid = np.linspace(x0, x1, nbins + 1)
        zpcorr = {}
        for band in self.bands:
            corrfile = os.path.join(self._path,
                        f"assets/zpcorr_idr3/SPLUS_{band}_offsets_grid.npy")
            corr = np.load(corrfile)
            zpcorr[band] = RectBivariateSpline(xgrid, xgrid, corr)
        return zpcorr

    def make_cube(self):
        fields = [_.replace("_", "-") for _ in self.zps["FIELD"]]
        flams, flamerrs = [], []
        headers = []
        for band, img in zip(self.bands, self.cutnames):
            wave = self.wave_eff[band] * u.Angstrom
            filename = os.path.join(self.cutouts_dir, img)
            h = fits.getheader(filename, ext=1)
            tile = h["OBJECT"].replace("_", "-")
            # Getting zero-point calibration and location correction
            zp0 = self.zps[fields.index(tile)][band]
            x0 = h["X0TILE"]
            y0 = h["Y0TILE"]
            zpcorr = float(self.zpcorr[band](x0, y0))
            zp = zp0 + zpcorr
            gain = h["GAIN"]
            f0 = np.power(10, -0.4 * (48.6 + zp))
            # Calculating flux density
            data = fits.getdata(filename, 1)
            fnu = data * f0 * self.fnu_unit
            flam = fnu * const.c / wave**2
            flam = flam.to(self.flam_unit).value
            # Uncertaintis in flux density
            weights = fits.getdata(filename.replace(".fz", "_weight.fz"), 1)
            dataerr = 1 / weights + np.clip(data, 0, np.infty) / gain
            fnuerr = dataerr * f0 * self.fnu_unit
            flamerr = fnuerr * const.c / wave**2
            flamerr = flamerr.to(self.flam_unit).value
            flams.append(flam)
            flamerrs.append(flamerr)
            headers.append(h)
        flam = np.array(flams)
        flamerr = np.array(flamerrs)
        # Making new header with WCS
        wcs = WCS(h)
        add_before_ind = 2
        inds = [i + 1 for i in range(wcs.wcs.naxis)]
        inds.insert(add_before_ind, 0)
        newwcs = wcs.sub(inds)
        newwcs.wcs.ctype[add_before_ind] = ''
        newwcs.wcs.cname[add_before_ind] = ''
        # Making new header template
        newheader = headers[0].copy()
        newheader.update(newwcs.to_header())
        # Making table with metadata
        tab = []
        tab.append(self.bands)
        tab.append([self.wave_eff[band] for band in self.bands])
        names = ["FILTER", "WAVE_EFF"]
        hfields = ["GAIN", "PSFFWHM", "DATE-OBS", "EXPTIME",
                   "EFECTIME", "NCOMBINE", "HIERARCH OAJ PRO FWHMMEAN"]
        for f in hfields:
            if not all([f in h for h in headers]):
                continue
            tab.append([h[f] for h in headers])
            names.append(f)
            if f in newheader:
                del newheader[f]
        tab = Table(tab, names=names)
        tab.rename_column("HIERARCH OAJ PRO FWHMMEAN", "PSFFWHM")

        # Producing data cubes HDUs.
        hdus = [fits.PrimaryHDU()]
        hdu1 = fits.ImageHDU(flam, newheader)
        hdu1.header["EXTNAME"] = ("DATA", "Name of the extension")
        hdus.append(hdu1)
        hdu2 = fits.ImageHDU(flamerr, newheader)
        hdu2.header["EXTNAME"] = ("ERRORS", "Name of the extension")
        hdus.append(hdu2)
        for hdu in hdus:
            hdu.header["BSCALE"] = (1, "Linear factor in scaling equation")
            hdu.header["BZERO"] = (0, "Zero point in scaling equation")
            hdu.header["BUNIT"] = ("{}".format(self.flam_unit),
                                   "Physical units of the array values")
        thdu = fits.BinTableHDU(tab)
        hdus.append(thdu)
        thdu.header["EXTNAME"] = "METADATA"
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(self.cubename, overwrite=True)
