#!/usr/bin/env python

# Python2 compatibility
from __future__ import print_function, division, absolute_import

import os
import sys
import numpy as np
from astropy.io import fits as astrofits
import argparse
import pytz
from datetime import datetime, timedelta
import logging
from tkp.accessors.dataaccessor import DataAccessor
from tkp.utility.coordinates import WCS

from lsl.common.mcs import mjdmpm_to_datetime

from lsl_toolkits.PasiImage import PasiImageDB

logger = logging.getLogger(__name__)

class pimsImage(DataAccessor):
    """Use the LWA Software Library to pull image data out of a pims file.
    Provide standard attributes, as per :class:`DataAccessor`."""
    def __init__(self, url, time_index, plane=None, beam=None):
        super(pimsImage, self).__init__()
        self.url = url
        try:
            db = PasiImageDB(self.url, mode='r')
        except Exception as e:
            print("ERROR: %s" % str(e))
        self.header, self.dbdata, self.spec = db.readImage() 
        self.data = self.read_data(time_index,self.dbdata)
        self.imSize = data.shape[-1]
        pScale = self.header['xPixelSize']
        self.sRad   = 360.0/pScale/np.pi / 2
        self.wcs = self.parse_coordinates()
        self.taustart_ts, self.tau_time = self.parse_times()
        self.freq_eff, self.freq_bw = self.parse_frequency()
        self.pixelsize = self.parse_pixelsize()
        bmaj,bmin,bpa = self.parse_beam(self.freq_eff)
        self.beam = self.degrees2pixels(
            bmaj, bmin, bpa, self.pixelsize[0], self.pixelsize[1]
            )
        self.centre_ra, self.centre_decl = self.calculate_phase_centre
        db.close()
        

    def _get_header(self):
        header, data, spec = self.db.readImage() 
        return header
    
    def read_data(self, time_index, dbdata)
        currentdata = dbdata[time_index]
        n_dim = len(data.shape)
        if plane is not None:
            data = data[plane]
        elif n_dim != 2:
            logger.warning("Loaded datacube with %s dimensions, assuming Stokes I and taking plane 0" % n_dim)
            data = data[0, :, :]
        data = data.transpose()

    def parse_coordinates(self):
        """Returns a WCS object"""
        header = self.header
        sRad = self.sRad
        imSize = self.imSize
        wcs = WCS()
        try:
            
            wcs.crval = header['zenithRA'],header['zenithDec']
            wcs.crpix = imSize/2 + 1 + 0.5 * ((imSize+1)%2),imSize/2 + 1 + 0.5 * ((imSize+1)%2)
            wcs.cdelt = -360.0/(2*sRad)/numpy.pi, 360.0/(2*sRad)/numpy.pi
        except KeyError:
            msg = "Coordinate system not specified in pims"
            logger.error(msg)
            raise TypeError(msg)
        wcs.ctype = 'RA---SIN','DEC--SIN'
        wcs.crota = 0., 0.
        wcs.cunit = header['cunit1'], header['cunit2']
        wcs.cunit = 'deg', 'deg'
        return wcs

    def calculate_phase_centre(self):
        x, y = self.data.shape
        centre_ra, centre_decl = self.wcs.p2s((x / 2, y / 2))
        return float(centre_ra), float(centre_decl)

    def parse_frequency(self):
        """
        Set some 'shortcut' variables for access to the frequency parameters
        in the FITS file header.

        @param hdulist: hdulist to parse
        @type hdulist: hdulist
        """
        freq_eff = self.header['freq']
        freq_bw = self.header['bandwidth']
        return freq_eff, freq_bw

    def parse_beam(self):
        """Read and return the beam properties bmaj, bmin and bpa values from
        the fits header.

        Returns:
          - Beam parameters, (semimajor, semiminor, position angle)
            in (pixels, pixels, radians)
        """
        header = self.header
        beamSize = 2.2*74e6/header['freq'] # lwa1
        bmaj, bmin, bpa = beamSize/header['xPixelSize'], beamSize/header['xPixelSize'],0.0


        return bmaj, bmin, bpa


    def parse_start_time(self):
        """
        Returns:
          - start time of image as an instance of ``datetime.datetime``
        """
        header = self.header
        mjd = int(header['startTime'])
        mpm = int((header['startTime'] - mjd)*86400.0*1000.0)
        tInt = header['intLen']*86400.0
        start = mjdmpm_to_datetime(mjd, mpm)
        return start



    def parse_times(self):
        """Returns:
          - taustart_ts: tz naive (implicit UTC) datetime at start of observation.
          - tau_time: Integration time, in seconds
        """
        # Attempt to do something sane with timestamps.
        header= self.header
        tInt = header['intLen']*86400.
        start = self.parse_start_time()
        end = start + timedelta(seconds=int(tInt), microseconds=int((tInt-int(tInt))*1000000))

        delta = end - start
        tau_time = delta.total_seconds()

        #For simplicity, the database requires naive datetimes (implicit UTC)
        #So we convert to UTC and then drop the timezone:
        timezone = pytz.timezone('US/Mountain') # LWA is in the Mountain timezone
        start_w_tz = start.replace(tzinfo=timezone)
        start_utc = pytz.utc.normalize(start_w_tz.astimezone(pytz.utc))
        return start_utc.replace(tzinfo=None), tau_time

      
