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

from tkp.accessors.OrvilleImageDB import OrvilleImageDB

logger = logging.getLogger(__name__)

class oimsImage(DataAccessor):
    """Use the LWA Software Library to pull image data out of a oims file.
    Provide standard attributes, as per :class:`DataAccessor`."""
    def __init__(self, filename, time_index, beamfile, freqind,plane=None):
        super(oimsImage, self).__init__()
        self.url = f"{filename},{time_index},{beamfile},{freqind}"
        try:
            db = OrvilleImageDB(filename,'r')
        except Exception as e:
            print("ERROR: %s" % str(e))
        db.seek(time_index)
        self.header, self.alldata = db.read_image()
        self.freqind = freqind
#         self.pb = np.load(beamfile)
        self.data = np.transpose(self.read_data(plane))# /self.pb)
        self.imSize = self.data.shape[-1]
        pScale = self.header['pixel_size']
        self.sRad  = 360.0/pScale/np.pi / 2
        self.wcs = self.parse_coordinates()
        self.taustart_ts, self.tau_time = self.parse_times()
        self.freq_eff, self.freq_bw = self.parse_frequency()
        self.pixelsize = self.parse_pixelsize()
        bmaj,bmin,bpa = self.parse_beam()
        self.beam = self.degrees2pixels(
            bmaj, bmin, bpa, self.pixelsize[0], self.pixelsize[1]
            )
        self.centre_ra, self.centre_decl = self.calculate_phase_centre()
        self.telescope="LWA-SV"
        db.close()
        
   
    
    def read_data(self,plane):
        start = datetime.now()
        data = self.alldata[self.freqind]
        data = data[0, :, :]
        t1 = datetime.now()
        return data

    def parse_coordinates(self):
        """Returns a WCS object"""
        header = self.header
        sRad = self.sRad
        imSize = self.imSize
        wcs = WCS()
        try:
            
            wcs.crval = header['center_ra'],header['center_dec']
            wcs.crpix = imSize/2 + 1 + 0.5 * ((imSize+1)%2),imSize/2 + 1 + 0.5 * ((imSize+1)%2)
            wcs.cdelt = -360.0/(2*sRad)/np.pi, 360.0/(2*sRad)/np.pi
        except KeyError:
            msg = "Coordinate system not specified in pims"
            logger.error(msg)
            raise TypeError(msg)
        wcs.ctype = 'RA---SIN','DEC--SIN'
        wcs.crota = 0., 0.
        wcs.cunit = 'deg', 'deg'
        return wcs

    def calculate_phase_centre(self):
        return self.header['center_ra'], self.header['center_dec']

    def parse_frequency(self):
        """
        Set some 'shortcut' variables for access to the frequency parameters
        in the FITS file header.

        @param hdulist: hdulist to parse
        @type hdulist: hdulist
        """
        hdr = self.header
        midfreq = (hdr['start_freq']  + ((self.freqind+1)*hdr['bandwidth']/2))
        freq_eff = midfreq 
        freq_bw = hdr['bandwidth']
        return freq_eff, freq_bw

    def parse_beam(self):
        """Read and return the beam properties bmaj, bmin and bpa values from
        the fits header.

        Returns:
          - Beam parameters, (semimajor, semiminor, position angle)
            in (pixels, pixels, radians)
        """
        hdr = self.header
        psize = hdr['pixel_size']
        midfreq = (hdr['start_freq']  + ((self.freqind+1)*hdr['bandwidth']/2))
        beamSize = 2.2*74e6/midfreq # lwa1
        bmaj, bmin, bpa = beamSize/psize, beamSize/psize,0.0


        return bmaj, bmin, bpa


    def parse_start_time(self):
        """
        Returns:
          - start time of image as an instance of ``datetime.datetime``
        """
        hdr = self.header
        mjd = int(hdr['start_time'])
        mpm = int((hdr['start_time'] - mjd)*86400.0*1000.0)
        tInt = hdr['int_len']*86400.0
        start = mjdmpm_to_datetime(mjd, mpm)
        return start



    def parse_times(self):
        """Returns:
          - taustart_ts: tz naive (implicit UTC) datetime at start of observation.
          - tau_time: Integration time, in seconds
        """
        # Attempt to do something sane with timestamps.
        hdr = self.header
        mjd = int(hdr['start_time'])
        mpm = int((hdr['start_time'] - mjd)*86400.0*1000.0)
        tInt = hdr['int_len']*86400.0
        start = self.parse_start_time()
        end = start + timedelta(seconds=int(tInt), microseconds=int((tInt-int(tInt))*1000000))

        delta = end - start
        tau_time = delta.total_seconds()

        #For simplicity, the database requires naive datetimes (implicit UTC)
        #So we convert to UTC and then drop the timezone:
        # oims images are in UTC
        return start, tau_time

      
