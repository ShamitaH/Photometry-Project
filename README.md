# Photometry Project

## Introduction
This Python project automates the processing of astronomical FITS files, including calibration, correction, and analysis, using the libraries `astropy`, `ccdproc`, and `photutils`. Additionally, it integrates with the Astrometry.net API to align and calibrate images based on celestial coordinates.

## Features
- **Image Calibration**: Automatic combination of multiple flat and dark frames to refine image quality.
- **Astrometry Correction**: Integration with Astrometry.net for accurate celestial coordinate calibration and WCS (World Coordinate System) information.
- **Star Detection**: Detection and analysis of stars using `DAOStarFinder`.
- **Aperture Photometry**: Measurement of star brightness through aperture photometry techniques.
- **Visual Output**: Generation of detailed visualizations showing processed images with highlighted detected stars.

## Dependencies
- Python 3.x
- numpy
- astropy
- ccdproc
- photutils
- matplotlib
- requests
- json (included in the Python standard library)
- time (included in the Python standard library)

Install the required Python packages using:
```bash
pip install numpy astropy ccdproc photutils matplotlib requests
