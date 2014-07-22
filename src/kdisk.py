# kdisk.py --- Compute a synthetic datacube for a Keplerian disk

import argparse
import configparser
from numpy import *
from astropy.constants import G, M_sun, au
from astropy.io import fits

def main():

    # Read options and arguments

    parser = argparse.ArgumentParser(description = "Compute a synthetic datacube for a Keplerian disk")
    parser.add_argument("parfile", help = "parameter file")
    args = parser.parse_args()
    
    # Read the parameters file

    params = configparser.ConfigParser()
    params.read(args.parfile)

    mstar = params.getfloat("disk", "mstar")
    size = params.getfloat("disk", "size")
    incl = params.getfloat("disk", "incl")
    pa = params.getfloat("disk", "pa")
    dist = params.getint("disk", "dist")
    lineint = params.get("line", "intensity").split(",")
    linewidth = params.getfloat("line", "width") 
    frequency = params.getfloat("line", "frequency")
    vlsr = params.getfloat("line", "vlsr")
    npix = params.getint("cube", "npix")
    pixsize = params.getfloat("cube", "pixsize")
    nchan = params.getint("cube", "nchan")
    chanwidth = params.getfloat("cube", "chanwidth")
    fitsname = params.get("output", "name")

    # Compute the initial grid

    ra_offset = (arange(npix) - (npix - 1) / 2) * pixsize
    dec_offset = ra_offset
    ra_grid, dec_grid = meshgrid(ra_offset, dec_offset)
    veloc = (arange(nchan) - (nchan - 1) / 2) * chanwidth + vlsr

    # Convert the projected coordinates to cylindrical coordinates in the plane of the disk

    pa *= pi / 180. # degrees -> rad
    x_grid = (-ra_grid * cos(pa - pi/2) - dec_grid * sin(pa - pi/2)) # disk major axis
    y_grid = (-ra_grid * sin(pa - pi/2) + dec_grid * cos(pa - pi/2)) # disk minor axis
    theta = zeros((npix, npix)) # theta = 0 along the l.o.s.
    r = zeros((npix, npix))
    mask = y_grid != 0
    incl *= pi / 180. # degrees -> rad
    theta[mask] = 2 * arctan((y_grid[mask] / cos(incl)) \
                             / (x_grid[mask] \
                                + sqrt(x_grid[mask]**2 + (y_grid[mask] / cos(incl))**2)))
    r = sqrt(x_grid**2 + (y_grid / cos(incl))**2) * dist # AU

    # Compute the line peak intensity

    peakint = zeros((npix, npix))
    if lineint[0] == "powerlaw":
        int_r1, r1, int_expn = map(lambda x: float(x), lineint[1:4]) 
        mask = r != 0
        peakint[mask] = int_r1 * pow(r[mask] / r1, int_expn)
    else:
        int0, fwhm = map(lambda x: float(x), lineint[1:3])
        sigma = fwhm / (2 * sqrt(2 * log(2)))
        peakint = int0 * exp(-r / (2 * sigma**2) / dist)
    if size != 0:
        peakint[r > size] = 0.

    # Compute the projection of the Keplerian velocity along the line of sight

    vproj = zeros((npix, npix))
    mask = r != 0
    vproj[mask] = sin(incl) * cos(theta[mask]) * sqrt(G * mstar * M_sun / (r[mask] * au))
    vproj *= 1e-3 # m/s -> km/s
    vproj[mask] += vlsr

    # Compute the synthetic datacube
    
    intensity = zeros((nchan, npix, npix))
    sigma = linewidth / (2 * sqrt(2 * log(2)))
    intensity = peakint[newaxis,:,:] *  exp(-(veloc[:,newaxis,newaxis] - vproj)**2 / (2 * sigma**2))

    # Export the datacube to FITS
    
    hdu = fits.PrimaryHDU(intensity)
    hdr = hdu.header
    hdr["CTYPE1"] = "RA---SIN"
    hdr["CTYPE2"] = "DEC--SIN"
    hdr["CTYPE3"] = "VELO-LSR"
    hdr["CUNIT1"] = "DEG"
    hdr["CUNIT2"] = "DEG"
    hdr["CUNIT3"] = "M/S"
    hdr["CDELT1"] = -pixsize / 3600.
    hdr["CDELT2"] = pixsize / 3600.
    hdr["CDELT3"] = chanwidth * 1e3
    hdr["CRPIX1"] = (npix - 1) / 2. + 1  # Fortran indexing
    hdr["CRPIX2"] = (npix - 1) / 2. + 1
    hdr["CRPIX3"] = (nchan - 1) / 2. + 1
    hdr["CRVAL1"] = 0.
    hdr["CRVAL2"] = 0.
    hdr["CRVAL3"] = vlsr * 1e3
    hdr["BUNIT"] = "K"
    hdr["RESTFREQ"] = frequency*1e6
    hdr["VELO-LSR"] = vlsr * 1e3
    hdu.writeto("%s.fits" % fitsname, clobber = True)

if __name__ == '__main__':
    main()
