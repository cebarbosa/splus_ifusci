""" Tool to create RGB images. """
import os

import numpy as np
from astropy.stats import sigma_clipped_stats
from PIL import Image
import cv2

def make_RGB_with_overlay(r, g, b, output, overlay=None, bb=5):
    """ Tool to produce RGB images with overlay.

    Parameters
    ----------
    r, g, b: numpy.array
        Arrays with the fluxes to be used in the R, G and B channels of the
        image.
    output: str
        Name of the output file.
    overlay: numpy.array
        Image to be overlay over RGB image.
    bb: float
        Smoothness parameter.
    """
    r[~np.isfinite(r)] = 0.
    g[~np.isfinite(g)] = 0.
    b[~np.isfinite(b)] = 0.
    I = (b + g + r) / 3.
    beta = np.nanmedian(I) * bb
    R = r * np.arcsinh(I / beta) / I
    G = g * np.arcsinh(I / beta) / I
    B = b * np.arcsinh(I / beta) / I
    maxRGB = np.percentile(np.stack([R, G, B]), 99.5)
    R = np.clip(255 * R / maxRGB, 0, 255).astype("uint8")
    G = np.clip(255 * G / maxRGB, 0., 255).astype("uint8")
    B = np.clip(255 * B / maxRGB, 0, 255).astype("uint8")
    RGB = np.stack([np.rot90(R, 3), np.rot90(G, 3),
                    np.rot90(B, 3)]).T
    if overlay is None:
        out = Image.fromarray(RGB)
        out.save(output)
        return
    mean, median, stddev = sigma_clipped_stats(overlay)
    maxha = np.percentile(overlay, 99.5)
    overlay = np.clip(overlay, stddev, maxha)
    overlay -= overlay.min()
    overlay = (overlay / overlay.max() * 255).astype("uint8")
    hamask = np.zeros_like(RGB)
    hamask[:, :, 0] = np.rot90(overlay, 3).T
    out = Image.fromarray(cv2.addWeighted(RGB, 0.9, hamask, 0.4, 0))
    out.save(output)