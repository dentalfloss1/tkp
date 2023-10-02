"""
This module provides the PasiImageDB class, which manages transactions
with a binary file format that stores PASI images.
"""

# Python2 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    range = xrange
    
import os
import sys
import ctypes
import shutil
import struct
import numpy as np

__version__ = '0.2'
__all__ = ['PasiImageDB',]


class _PrintableLittleEndianStructure(ctypes.LittleEndianStructure):
    """
    Sub-class of ctypes.LittleEndianStructure that adds a as_dict()
    method for accessing all of the fields as a dictionary.
    """
    
    def __getitem__(self, key):
        value = getattr(self, key, None)
        if value is None:
            raise KeyError("Cannot find key '%s'" % key)
        return value
        
    def __setitem__(self, key, value):
        setattr(self, key, value)
        
    def keys(self):
        return [field[0] for field in self._fields_]
        
    def __contains__(self, value):
        return value in self.keys()
        
    def sizeof(self):
        return ctypes.sizeof(self)
        
    def as_dict(self):
        """
        Return all of the structure fields as a dictionary.
        """
        
        out = {}
        for field in self._fields_:
            out[field[0]] = getattr(self, field[0], None)
        return out
        
    def __repr__(self):
        return repr(self.as_dict())

class _PasiImageDBv001_FileHeader(_PrintableLittleEndianStructure):
    _pack_   = 1
    _fields_ = [('corrVersion', ctypes.c_char*16),
                ('imagerVersion', ctypes.c_char*16),
                ('station', ctypes.c_char*16),
                ('stokesParams', ctypes.c_char*16),
                ('xSize', ctypes.c_uint),
                ('ySize', ctypes.c_uint),
                ('nSpecChans', ctypes.c_uint),
                ('flags', ctypes.c_uint),
                ('startTime', ctypes.c_double),
                ('stopTime', ctypes.c_double)]

class _PasiImageDBv002_FileHeader(_PrintableLittleEndianStructure):
    _pack_   = 1
    _fields_ = [('corrVersion', ctypes.c_char*16),
                ('imagerVersion', ctypes.c_char*16),
                ('station', ctypes.c_char*16),
                ('stokesParams', ctypes.c_char*16),
                ('xSize', ctypes.c_uint),
                ('ySize', ctypes.c_uint),
                ('nSpecChans', ctypes.c_uint),
                ('flags', ctypes.c_uint),
                ('xPixelSize', ctypes.c_double),
                ('yPixelSize', ctypes.c_double),
                ('startTime', ctypes.c_double),
                ('stopTime', ctypes.c_double)]

_PasiImageDBv003_FileHeader = _PasiImageDBv002_FileHeader

class _PasiImageDBv001_IntHeader(_PrintableLittleEndianStructure):
    _pack_   = 1
    _fields_ = [('visFileName', ctypes.c_char*256),
                ('startTime', ctypes.c_double),
                ('centroidTime', ctypes.c_double),
                ('intLen', ctypes.c_double),
                ('lst', ctypes.c_double),
                ('freq', ctypes.c_double),
                ('bandwidth', ctypes.c_double),
                ('gain', ctypes.c_double),
                ('zenithRA', ctypes.c_double),
                ('zenithDec', ctypes.c_double),
                ('worldreplace0', ctypes.c_double*2)]

_PasiImageDBv002_IntHeader = _PasiImageDBv001_IntHeader

class _PasiImageDBv003_IntHeader(_PrintableLittleEndianStructure):
    _pack_   = 1
    _fields_ = [('visFileName', ctypes.c_char*256),
                ('startTime', ctypes.c_double),
                ('centroidTime', ctypes.c_double),
                ('intLen', ctypes.c_double),
                ('lst', ctypes.c_double),
                ('freq', ctypes.c_double),
                ('bandwidth', ctypes.c_double),
                ('gain', ctypes.c_double),
                ('fill', ctypes.c_double),
                ('zenithRA', ctypes.c_double),
                ('zenithDec', ctypes.c_double),
                ('worldreplace0', ctypes.c_double*2)]


class PasiImageDB(object):
    """
    Encapsulates a PasiImageDB binary file.
    
    This class can be used for both reading and writing PasiImageDB files.
    For reading, initialize with mode = "r" and use the readImage() function.  
    For writing, use mode = "a" or "w" and the addImage() function. Be sure 
    to always call the close() method after adding images in order to update 
    the file's header information.
    
    Public module variables:
      header -- a class with member variables describing the data, including:
        corrVersion, a string giving the correlator version
        imagerVersion, ditto for the imager
        station, the station name
        stokesParams, the comma-delimited parameters (e.g., 'I,Q,U,V')
        xSize and ySize, the dimensions of the image in pixels
        xPixelSize and yPixelSize, the physical size of a pixel, in degrees
        nSpecChans, the number of channels in the spectrum
        flags, a bitfield: 0x1 = sorted; others zero
        startTime, the earliest time of data covered by this file, in MJD UTC
        stopTime, the latest time of data, in MJD UTC
      file -- the underlying file object, which probably shouldn't be touched
      version -- the format version of the output file
    """
    
    # The PasiImageDB files start with a 16 byte string that specifies the
    # data file's format version.  Following it is that format's header block,
    # defined by the _fileHeaderStructs dictionary.  After the header are the
    # integration blocks.  An integration block starts with a header, defined
    # by the _intHeaderStructs dictionary.  Following that is a single
    # dipole's spectrum, encoded in float32s, and then the images, packed as
    # [stokes, x, y] and also float32s.  The alignment of the spectrum and
    # images is machine-dependent due to me being lazy and using
    # numpy.fromfile() and ndarray.tofile(); all others are low-endian.
    #
    # All absolute times are in MJD UTC, LSTs are in days (i.e., from 0 to
    # 0.99726957), and integration lengths are in days.  Sky directions
    # (including RA) and pixel sizes are in degrees.  All other entries are in
    # standard mks units.
    
    _currentFormatVersion = b'PasiImageDBv003'
    
    _fileHeaderStructs = {
        b'PasiImageDBv001': _PasiImageDBv001_FileHeader,
        b'PasiImageDBv002': _PasiImageDBv002_FileHeader,
        b'PasiImageDBv003': _PasiImageDBv003_FileHeader,
        }
    flagSorted = 0x0001
    
    _intHeaderStructs = {
        b'PasiImageDBv001': _PasiImageDBv001_IntHeader,
        b'PasiImageDBv002': _PasiImageDBv002_IntHeader,
        b'PasiImageDBv003': _PasiImageDBv003_IntHeader
        }
    _timeOffsets = {b'PasiImageDBv001': 256,
                    b'PasiImageDBv002': 256,
                    b'PasiImageDBv003': 256 }
    
    def __init__(self, fileName, mode = 'r',
                 corrVersion = '', imagerVersion = '', station = ''):
        """
        Constructs a new PasiImageDB.
        
        Optional arguments specify the file mode (must be 'r', 'w', or, 'a';
        defaults to 'r') and strings providing the correlator version, the
        imager version, and the station name, all of which are truncated at 16
        bytes.  These optional strings are only relevant when opening a file
        for writing.
        """
        
        self.file = None
        self.iIntegration = -1
        
        self._FileHeader = self._fileHeaderStructs[self._currentFormatVersion]
        self._IntHeader = self._intHeaderStructs[self._currentFormatVersion]
        self._TimeOffsets = self._timeOffsets[self._currentFormatVersion]
        
        # For read mode, we do not create a new file.  Raise an error if it
        # does not exist, and create an empty PasiImageDB object if its length
        # is zero.
        if mode == 'r':
            self._isNewFile = False
            if not os.path.isfile(fileName):
                raise OSError('The specified file, "%s", does not exist.'
                            % fileName)
            fileSize = os.path.getsize(fileName)
            if fileSize == 0:
                self.version = self._currentFormatVersion
                self.header = self._FileHeader()
                self.iIntegration = 0
                self.nIntegrations = 0
                self.nStokes = 0
                return
        
        # For append mode, check if the file exists and is at least longer
        # than the initial 16 byte version string.  If that's the case, switch
        # to 'r+' mode, since we may need to read and/or write to the header,
        # and some Unix implementations don't allow this with 'a' mode.
        # Otherwise, switch to write mode.
        elif mode == 'a':
            fileSize = os.path.getsize(fileName) \
                if os.path.isfile(fileName) else 0
            self._isNewFile = (fileSize <= 16)
            mode = 'w' if self._isNewFile else 'r+'
        
        # Write mode: pretty straightforward.
        elif mode == 'w':
            self._isNewFile = True
        
        else:
            raise ValueError("Mode must be 'r', 'w', or 'a'.")
        
        # Now read or create the file header.
        mode += 'b'
        self.file = open(fileName, mode)
        self._fileHeaderOutdated = False
        
        if not self._isNewFile:
            # For existing files, get the number of iterations from the file
            # length, and read the format version and file header.
            self.version = self.file.read(16).rstrip(b'\x00')
            if self.version not in self._fileHeaderStructs:
                raise KeyError('The file "%s" does not appear to be a '
                            'PasiImageDB file.  Initial string: "%s"' %
                            (fileName, self.version))
            self._FileHeader = self._fileHeaderStructs[self.version]
            self._IntHeader = self._intHeaderStructs[self.version]
            self._TimeOffsets = self._timeOffsets[self.version]
            headerStruct = self._FileHeader()
            
            if mode != 'r' and fileSize <= 16 + headerStruct.sizeof():
                # If the file is too short to have any data in it, close it
                # and start a new one.  This one is probably corrupt.
                self.file.close()
                self._isNewFile = True
                mode = 'w'
                self.file = open(fileName, mode)
            
            else:
                # It looks like we should have a good header, at least ....
                self.header = self._FileHeader()
                self.file.readinto(self.header)
                if self.version < b'PasiImageDBv002':
                    self.header.xPixelSize = 1.0  # Default to 1 deg/pix.
                    self.header.yPixelSize = 1.0
                self.nStokes = len(self.header.stokesParams.split(b','))
                
                intHeader = self._IntHeader()
                intSize = intHeader.sizeof() + \
                    4 * (self.header.nSpecChans +
                        self.nStokes * self.header.xSize * self.header.ySize)
                if (fileSize - 16 - self.header.sizeof()) % intSize != 0:
                    raise RuntimeError('The file "%s" appears to be '
                                    'corrupted.' % fileName)
                self.nIntegrations = \
                    (fileSize - 16 - self.header.sizeof()) // intSize
                
                if mode == 'r+b':
                    self.file.seek(0, os.SEEK_END)
                    self.iIntegration = self.nIntegrations
                else:
                    self.iIntegration = 0
        
        if self._isNewFile:
            # Start preparing a file header, but don't write it until we
            # receive the first image, which will fill in some information
            # (e.g., resolution) that isn't yet available.
            self.version = self._currentFormatVersion
            self.header = self._FileHeader()
            self.header.corrVersion = corrVersion
            self.header.imagerVersion = imagerVersion
            self.header.station = station
            self.header.flags = self.flagSorted
            self.nIntegrations = 0
    
    
    def __del__(self):
        if self.file is not None and not self.file.closed:
            self.close()
    
    
    def close(self):
        """
        Closes the database file.  If the header information is outdated, it
        writes the new file header.
        """
        
        if self.file is None or self.file.closed:  return
        
        if self._fileHeaderOutdated:
            self.file.seek(16, os.SEEK_SET)
            self.file.write(self.header)
        
        self.file.close()
        self.iIntegration = -1
    
    
    def closed(self):
        return self.file is None or self.file.closed
    
    
    def getpos(self):
        return self.iIntegration
    
    
    def eof(self):
        return self.iIntegration >= self.nIntegrations
    
    
    def seek(self, index):
        if index < 0:
            index += self.nIntegrations
        if index < 0 or index >= self.nIntegrations:
            raise IndexError('PasiImageDB index %d outside of range [0, %d)' %
                            (index, self.nIntegrations))
        if self.iIntegration != index:
            intHeader = self._IntHeader()
            intSize = intHeader.sizeof() + \
                4 * (self.header.nSpecChans +
                    self.nStokes * self.header.xSize * self.header.ySize)
            fileHeader = self._FileHeader()
            headerSize = 16 + fileHeader.sizeof()
            self.file.seek(headerSize + intSize * index, os.SEEK_SET)
            self.iIntegration = index
    
    
    def _checkHeader(self, stokesParams, xSize, ySize,
                    xPixelSize, yPixelSize, nSpecChans,
                    station = None):
        """
        For new files, adds the given information to the file header and
        writes the header to disk.  For existing files, compares the given
        information to the expected values and raises a ValueError if there's
        a mismatch.
        """
        
        if type(stokesParams) is list:
            stokesParams = ','.join(stokesParams)
        
        if self._isNewFile:
            # If this is the file's first image, fill in values of the file
            # header based on the image properties, then write the header.
            if station:
                self.header.station  = station
            self.header.stokesParams = stokesParams
            self.header.xSize        = xSize
            self.header.ySize        = ySize
            self.header.xPixelSize   = xPixelSize
            self.header.yPixelSize   = yPixelSize
            self.header.nSpecChans   = nSpecChans
            self.file.write(struct.pack('16s', self.version))
            headerStruct = self._fileHeaderStructs[self.version]()
            self.file.write(headerStruct)
            self.nStokes = len(self.header.stokesParams.split(b','))
            self._isNewFile = False
        
        else:
            # Make sure the station name (if given) matches expectations.
            if station and station != self.header.station:
                raise ValueError(
                    'The station given for this image ("%s") does not match '
                    'this file\'s station ("%s").' %
                    (station, self.header.station))
            
            # Make sure that the Stokes parameters match expectations.
            if stokesParams != self.header.stokesParams:
                raise ValueError(
                    'The Stokes parameters for this image (%s) do not match '
                    'this file\'s parameters (%s).' %
                    (stokesParams, self.header.stokesParams))
            
            # Make sure that the dimensions of the data match expectations.
            if xSize != self.header.xSize or ySize != self.header.ySize:
                raise ValueError(
                    'The spatial resolution of this image (%d x %d) does not '
                    'match this file\'s resolution (%d x %d).' %
                    (xSize, ySize, self.header.xSize, self.header.ySize))
            
            if xPixelSize != self.header.xPixelSize or \
            yPixelSize != self.header.yPixelSize:
                raise ValueError(
                    'The pixel size of this image (%r deg x %r deg) does not '
                    'match this file\'s resolution (%r deg x %r deg).' %
                    (xPixelSize, yPixelSize,
                    self.header.xPixelSize, self.header.yPixelSize))
            
            # Make sure that the size of the spectrum matches expectations.
            if nSpecChans != self.header.nSpecChans:
                raise ValueError(
                    'The length of the spectrum for this image (%d) does not '
                    'match this file\'s spectrum length (%d).'
                    % (nSpecChans, self.header.nSpecChans))
    
    
    def _updateHeaderInfo(self, interval):
        """
        To be called at the end of the addImage functions.  Updates the header
        information to reflect the new data.
        """
        self.nIntegrations += 1
        
        # Has this image expanded the time range covered by the file?
        if self.header.startTime == 0 or \
        self.header.startTime > interval[0]:
            self.header.startTime = interval[0]
            self._fileHeaderOutdated = True
        
        if self.header.stopTime < interval[1]:
            self.header.stopTime = interval[1]
            self._fileHeaderOutdated = True
        
        # If the new image isn't later than all the others, and the file is
        # currently marked as sorted, then remove the sorted flag.
        elif self.header.flags & self.flagSorted:
            self.header.flags &= ~self.flagSorted
            self._fileHeaderOutdated = True
    
    def addImage(self, info, data, spec = None):
        """
        Adds an integration to the database.  Returns the index of the newly
        added image.
        
        Arguments:
          info -- an object with the following member variables defined:
            visFileName -- the name of the visibility file that was imaged
            startTime -- MJD UTC at which this integration began
            centroidTime -- mean MJD UTC of this integration
            intLen -- integration length, in days
            lst -- mean local sidereal time of the observation, in days
            freq -- center frequency of the integration, in Hz
            bandwidth -- bandwidth of the integrated data, in Hz
            gain -- gain setting of the observation (-1 if unknown)
            fill -- fill fraction of this image's integration (-1 if unknown)
            zenithRA -- RA of zenith, in degrees
            zenithDec -- Declination of zenith, in degrees
            xPixelSize -- Real-world width of a pixel, in degrees
            yPixelSize -- Real-world height of a pixel, in degrees
            worldreplace0 -- CASA-defined two-element float; purpose unknown
            stokesParams -- a list or comma-delimited string of Stokes params
          data -- a 3D float array of image data indexed as [iStokes, x, y]
          spec -- an optional array containing the single dipole spectrum
        """
        
        self._checkHeader(
            info.stokesParams, data.shape[1], data.shape[2],
            info.xPixelSize, info.yPixelSize,
            len(spec) if spec is not None else 0)
        
        # Write it out.
        intHeader = self._IntHeader()
        for key in info.keys():
            if key not in intHeader.keys():
                continue
            setattr(intHeader, key, info[key])
        self.file.write(intHeader)
        if spec is not None:
            spec.astype(np.float32).tofile(self.file)
        data.astype(np.float32).tofile(self.file)
        self.file.flush()
        interval = [info.startTime, info.startTime + info.intLen]
        self._updateHeaderInfo(interval)
        return self.nIntegrations - 1
    
    
    def readImage(self):
        """
        Reads an integration from the database.
        
        Returns a 3-tuple containing:
        info -- an object with the following member variables defined:
            visFileName -- the name of the visibility file that was imaged
            startTime -- MJD UTC at which this integration began
            centroidTime -- mean MJD UTC of this integration
            intLen -- integration length, in days
            lst -- mean local sidereal time of the observation, in days
            freq -- center frequency of the integration, in Hz
            bandwidth -- bandwidth of the integrated data, in Hz
            gain -- gain setting of the observation (-1 if unknown)
            fill -- fill fraction of the integration (-1 if unknown)
            zenithRA -- RA of zenith, in degrees
            zenithDec -- Declination of zenith, in degrees
            xPixelSize -- Real-world width of a pixel, in degrees
            yPixelSize -- Real-world height of a pixel, in degrees
            worldreplace0 -- CASA-defined two-element float; purpose unknown
            stokesParams -- a list or comma-delimited string of Stokes params
        data -- a 3D float array of image data indexed as [iStokes, x, y]
        spec -- if available, an array containing the single dipole spectrum
        """
        
        intHeader = self._IntHeader()
        self.file.readinto(intHeader)
        intHeader.xPixelSize = self.header.xPixelSize
        intHeader.yPixelSize = self.header.yPixelSize
        intHeader.stokesParams = self.header.stokesParams
        if 'gain' not in intHeader:
            intHeader.gain = -1
        if 'fill' not in intHeader:
            intHeader.fill = -1.
            
        nStokes, cx, cy = self.nStokes, self.header.xSize, self.header.ySize
        if self.header.nSpecChans > 0:
            spec = np.fromfile(self.file, np.float32, self.header.nSpecChans)
        else:
            spec = None
        data = np.fromfile(self.file, np.float32, nStokes * cx * cy
                        ).reshape(nStokes, cx, cy)
        
        self.iIntegration += 1
        return intHeader, data, spec
    
    
    @staticmethod
    def sort(fileName):
        """
        Sorts the integrations in a DB file to be time-ordered.
        """
        
        # Open the input database.  If it's already sorted, stop.
        inDB = PasiImageDB(fileName, 'r')
        if inDB.header.flags & PasiImageDB.flagSorted:
            inDB.close()
            return
        
        # Read the entire input database into memory.
        inIntHeaderStruct = PasiImageDB._intHeaderStructs[inDB.version]()
        headerSize = inIntHeaderStruct.sizeof()
        dataSize = 4 * (inDB.header.nSpecChans +
                        inDB.nStokes * inDB.header.xSize * inDB.header.ySize)
        intSize = headerSize + dataSize
        
        data = inDB.file.read()
        inDB.file.close()
        if len(data) != intSize * inDB.nIntegrations:
            raise RuntimeError('The file "%s" appears to be corrupted.' %
                            fileName)
        
        # Loop throught the input DB's images, saving their image times.
        # Determine the sort order of those times.
        times = np.array([
                struct.unpack_from('d', data, offset = i)[0] for i in
                range(PasiImageDB._timeOffsets[inDB.version],
                    intSize * inDB.nIntegrations, intSize)])
        
        intOrder = times.argsort()
        
        # Write the sorted file.  Note that we write it using the most recent
        # header version, which may differ from the version of the input file.
        # After writing the updated file header, loop through the intervals
        # and copy (in sorted order) from the input data to the output file.
        inDB.header.flags |= PasiImageDB.flagSorted
        
        outFile = open(fileName, 'w')
        outVersion = PasiImageDB._currentFormatVersion
        outFileHeaderStruct = PasiImageDB._fileHeaderStructs[outVersion]()
        outIntHeaderStruct = PasiImageDB._intHeaderStructs[outVersion]()
        outFile.write(struct.pack('16s', outVersion))
        outFile.write(inDB.header)
        
        for iOut in range(inDB.nIntegrations):
            i = intOrder[iOut] * intSize
            header = PasiImageDB._intHeaderStructs[outVersion]()
            ctypes.memmove(ctypes.addressof(header), data[i : i + headerSize], headerSize)
            outFile.write(header)
            outFile.write(data[i + headerSize : i + intSize])
        
        outFile.close()
    
    
    # Implement some built-ins to make reading images more "Pythonic" ...
    def __len__(self):
        return self.nIntegrations
        
        
    def __getitem__(self, index):
        self.seek(index)
        return self.readImage()
        
        
    def __iter__(self):
        return self
        
        
    def __next__(self):
        if self.iIntegration >= self.nIntegrations:
            raise StopIteration
        else:
            return self.readImage()
            
            
    def next(self):
        return self.__next__()