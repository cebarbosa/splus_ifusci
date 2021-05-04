import numpy as np
import astropy.units as u

class SCube():
    def __init__(self, obj, coords, size):
        self.obj = obj
        self.coors = coords
        self.size = size
        # General definitions
        self.ps = 0.55 * u.arcsec / u.pixel
        self.bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R',
                      'F660', 'I', 'F861', 'Z']
        self.bands_names = {'U': "$u$", 'F378': "$J378$", 'F395': "$J395$",
                            'F410': "$J410$", 'F430': "$J430$", 'G': "$g$",
                            'F515': "$J515$", 'R': "$r$", 'F660': "$J660$",
                            'I': "$i$", 'F861': "$J861$", 'Z': "$z$"}
        self.wave_eff = {"F378": 3770.0, "F395": 3940.0, "F410": 4094.0,
                         "F430": 4292.0, "F515": 5133.0, "F660": 6614.0,
                         "F861": 8611.0, "G": 4751.0, "I": 7690.0, "R": 6258.0,
                         "U": 3536.0, "Z": 8831.0}
        self.exptimes = {"F378": 660, "F395": 354, "F410": 177,
                         "F430": 171, "F515": 183, "F660": 870, "F861": 240,
                         "G": 99, "I": 138, "R": 120, "U": 681,
                         "Z": 168}
        self.wave = np.array([self.wave_eff[band] for band in \
                              self.bands]) *  u.Angstrom
        self.flam_unit = u.erg / u.cm / u.cm / u.s / u.AA
        self.fnu_unit = u.erg / u.s / u.cm / u.cm / u.Hz
        # Producing
