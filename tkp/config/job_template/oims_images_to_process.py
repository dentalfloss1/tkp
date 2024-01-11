
###################################################################################
#      List the images for processing by the transient detection pipeline         #
###################################################################################

# This should provide a module-scope iterable named "images" which provides a
# full path to each image to be processed. For example:

# images = [
#     "/path/to/image1",
#     "/path/to/image2",
# ]

# Optionally, whatever standard tools are required may be used to generate the
# list:
#
# import os
# import glob
# images = sorted(
#     glob.glob(
#         os.path.expanduser("/mnt/vg0/oims/*31.7*stokesI*fits")
#     )
# )

# If you want to read in pims files: uncomment the below lines
import os
import numpy as np 
import glob
from tkp.accessors.OrvilleImageDB import OrvilleImageDB
from lsl.common.mcs import mjdmpm_to_datetime
from lsl.common.paths import DATA as dataPath
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.wcs import WCS
from astropy.time import Time

def calcbeamprops(az,alt,header,freq):

    # az and alt need to be the same shape as the image we will correct

    i = 0
    beamDict = np.load(os.path.join(dataPath, 'lwa1-dipole-emp.npz'))
    polarpatterns = []
    for beamCoeff in (beamDict['fitX'], beamDict['fitY']):
        alphaE = np.polyval(beamCoeff[0,0,:],freq )
        betaE =  np.polyval(beamCoeff[0,1,:],freq )
        gammaE = np.polyval(beamCoeff[0,2,:],freq )
        deltaE = np.polyval(beamCoeff[0,3,:],freq )
        alphaH = np.polyval(beamCoeff[1,0,:],freq )
        betaH =  np.polyval(beamCoeff[1,1,:],freq )
        gammaH = np.polyval(beamCoeff[1,2,:],freq )
        deltaH = np.polyval(beamCoeff[1,3,:],freq )
        corrFnc = None

        def compute_beam_pattern(az, alt, corr=corrFnc):
            zaR = np.pi/2 - alt*np.pi / 180.0
            azR = az*np.pi / 180.0

            c = 1.0
            if corrFnc is not None:
                c = corrFnc(alt*np.pi / 180.0)
                c = np.where(np.isfinite(c), c, 1.0)

            pE = (1-(2*zaR/np.pi)**alphaE)*np.cos(zaR)**betaE + gammaE*(2*zaR/np.pi)*np.cos(zaR)**deltaE
            pH = (1-(2*zaR/np.pi)**alphaH)*np.cos(zaR)**betaH + gammaH*(2*zaR/np.pi)*np.cos(zaR)**deltaH

            return c*np.sqrt((pE*np.cos(azR))**2 + (pH*np.sin(azR))**2)
        # Calculate the beam
        pattern = compute_beam_pattern(az, alt)
        polarpatterns.append(pattern)
        i += 1
    beamDict.close()
    return polarpatterns[0], polarpatterns[1]

def pbcorr(pbfile,header,imSize,chan):
    mjd = int(header['start_time'])
    mpm = int((header['start_time'] - mjd)*86400.0*1000.0)
    tInt = header['int_len']*86400.0
    dateObs = mjdmpm_to_datetime(mjd, mpm)
    x = np.arange(imSize) - 0.5
    y = np.arange(imSize) - 0.5
    x,y = np.meshgrid(x,y)
    pScale = header['pixel_size']
    sRad   = 360.0/pScale/np.pi / 2
    crval1 = header['center_ra']*np.pi/180
    crpix1 = imSize/2 + 1 + 0.5 * ((imSize+1)%2) 
    cdelt1 = np.pi*(-360.0/(2*sRad)/np.pi)/180
    crval2 = header['center_dec']*np.pi/180
    crpix2 = imSize/2 + 1 + 0.5 * ((imSize+1)%2) 
    cdelt2 = np.pi*(360.0/(2*sRad)/np.pi)/180
    ra = ((crval1 + (x - crpix1)*cdelt1/(np.cos(crval2)))*180/np.pi) 
    dec = (crval2 + cdelt2*(y-crpix2))*180/np.pi
    # Make dec go between -90 and 90
    # Adjust RA accordingly
    decover = dec>90
    decdiff = dec[decover] - 90
    dec[decover] = dec[decover] - decdiff
    ra[decover] +=180
    decoverneg = dec<-90
    decdiffneg = dec[decoverneg] + 90
    dec[decoverneg] = dec[decoverneg] + decdiffneg
    ra[decoverneg] +=180
    ra = ra % 360
    
    sc = SkyCoord(ra,dec,unit='deg')
    lwasv = EarthLocation.from_geodetic(-106.885664,34.348562, height=1475) 
    time = Time(dateObs.strftime("%Y-%m-%dT%H:%M:%S"),format="isot")
    aa = AltAz(location=lwasv, obstime=time)
    myaltaz = sc.transform_to(aa)
    alt = myaltaz.alt.deg
    az = myaltaz.az.deg
    # Keep alt between 0 and 90, adjust az accordingly
    negalt = alt < 0
    alt[negalt] += 90
    az[negalt] + 180

    freq = (hdr['start_freq']  + ((chan+1)*hdr['bandwidth']/2))
    XX,YY = calcbeamprops(az,alt,header,freq)
    # Correct stokes I only, need XY and YX for U and V, could do Q here
    np.save(pbfile,((XX+YY)/2))

oimsfiles = glob.glob("/home/sarah/Data/oims/58956*.oims")
images = []
cwd = os.getcwd()
allfreq = []
for file in oimsfiles:
        db = OrvilleImageDB(file, mode='r') 
        ints = db.nint
        for i in range(ints):
            db.seek(i)
            hdr, alldata = db.read_image()
            for chan,stokesdata in enumerate(alldata):
                data = alldata[chan]
                freq = (hdr['start_freq']  + ((chan+1)*hdr['bandwidth']/2))
                pbpath = f"{cwd}/pb{freq}.npy"
                images.append(f"{file},{i},{pbpath},{chan}")
                allfreq.append((pbpath,hdr,data.shape[-1],chan))
        db.close()

allfreq = np.array(allfreq,dtype=[('file','<U256'),('header','O'),('imSize','i8'),('chan','i8')])
_,uniquefreqs = np.unique(allfreq['file'],return_index=True)
# for f in allfreq[uniquefreqs]:
#     pbcorr(*f)
#Display the list of images to be processed whenever this file is imported:
# (can be used for quick checking via an ipython import)
print("******** IMAGES: ********")
for f in images:
    print(f)
print("*************************")
