Kdisk
=====

`kdisk` is a simple Python program to compute the line emission from a
geometrically thin and optically thick protoplanetary disk. It creates
a datacube in FITS format that can be processed in a data reduction
package (such as GILDAS) to produce synthetic images and
visibilities. These synthetic data can be compared with observations
to determine the properties (e.g. central mass or inclination) of an
observed disk.

The disk is assumed to be in Keplerian rotation a radius lower than
the centrifugal radius (which can be set to a large value, for a
purely Keplerian disk), and in infall plus rotation beyond the
centrifugal radius.

Input file
----------

Model paramters are read from an input file. Here is an example:

```ini
[disk]
mstar = 0.20
incl = 85.
pa = 3.
r0 = 80
size = 0.
dist = 140
[cube]
npix = 1024
pixsize = 0.005
nchan = 128
chanwidth = 0.1
[line]
frequency = 219560.3541
intensity = gaussian, 135., 1.56
width = 0.1
vlsr = 6.0
[output]
name = l1527-c18o
```

Parameters
----------

- `mstar`:  mass of the central object in solar masses
- `incl`: disk inclination in degrees (90 for edge-on, 0 for face-on)
- `pa`: position angle of projected disk rotation axis (0 for a
  North-South, 90 for East-West)
- `rc`: centrifugal radius, in AU 
- `size`: disk size, in AU (set to 0 for an infinite disk)
- `dist`: disk distance, in pc
- `npix`: number of pixels in RA and Dec offset
- `pixsize`: size of a pixel, in arcsecs
- `nchan`: number of velocity channels
- `chanwidth`: channel width, in km/s
- `frequency`: line frequency, in MHz
- `intensity`: line intensity of the disk surface. For a Gaussian
  distribution, set to `gaussian,int0,fwhm` where `int0` is the peak
  intensity, and `fwmh` is the FWHM (in arcsecs). For a power-law, set
  to `powerlaw,int_r1,r1,int_expn` where `int_r1` is the intensity at
  the radius `r1` (in arcsecs), and `int_expn` is the powerlaw
  exponent. For a ring, set to `ring,int_ring,r1,r2`, where `int_ring`
  is the intensity between radii `r1` and `r2` (in arcsecs).
- `width`: line width, in km/s
- `vlsr`: source systemic velocity in the LSR, in km/s
- `name`: base name of the output FITS file

Usage
-----

```
% python kdisk.py input.ini`
```

where `input.ini` is the name of the output file.
