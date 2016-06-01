import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from numpy import zeros,linspace,where,size,arange,sum,min,int,ceil
from astropy import units as u
from astropy.time import Time
from datetime import datetime,timedelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.io import fits
from DARK_WINDOW import dark_window as dw

def GW_Followup(File_Name,WHEN='',min_obsv=30,nm=1,FIG=False,CLU=False):
    # comments:
    # min_obsv: should be always a factor of 10
    #
    tm_available = -10.	 
    night_strt, wh_sz, tm_available = dw(tmrw=0,WHEN=WHEN)
    if wh_sz <= min_obsv/10.: 
	night_strt, wh_sz, tm_available = dw(tmrw=1,WHEN=WHEN)
    #
    #################################################################
    File_Address = 'TEMP_'+File_Name+str(nm)+'sigma'+'.fits'	
    gwgc_slctd = fits.open(File_Address)
    gwgc_ra_1 = gwgc_slctd[1].data['RAJ2000']
    gwgc_dec_1 = gwgc_slctd[1].data['DEJ2000']
    gwgc_dec_allowed  = (gwgc_dec_1 >= -27.5) & (gwgc_dec_1 <= 57.)
    #
    gwgc_ra = gwgc_ra_1[gwgc_dec_allowed]
    gwgc_dec = gwgc_dec_1[gwgc_dec_allowed]
    #####
    if (CLU):
    	gwgc_M_B = gwgc_slctd[1].data['BMAG1']
	M_B = gwgc_M_B[gwgc_dec_allowed]
        ### CLU Galaxy Selection Criterion (L_B / L_B_star > x_(1/2) ) 
        LB_Lsun = 10**((M_B - 4.77)/-2.5).astype('float64')
        Lsun = 4.e33 # erg/s
        LB = LB_Lsun * Lsun
        LB_sun = 2.16e33 # erg/s
        #LB0 = 1.e-3 * 1.e10 * LB_sun
        h = 0.7
        LB_star = 1.2*1e10*(h**-2)*LB_sun
        x_cut = 0.626
        LB_CLU = x_cut * LB_star
        wh_CLU = (LB >= LB_CLU)
	gwgc_ra = gwgc_ra[wh_CLU]
	gwgc_dec = gwgc_dec[wh_CLU]
    #####
    #
    RA = gwgc_ra*u.deg
    DEC = gwgc_dec*u.deg
    ################################################################	    
    tw = 10 # minute -- time_window
    time_1 = night_strt
    delta_t = arange(0, wh_sz*tw, tw)*u.min
    time_1 = time_1 + delta_t
    
    k = 0
    strg = zeros((size(RA),wh_sz),dtype=bool)
    radecs = SkyCoord(RA,DEC)
    
    # Geodetic coordinates of observatory (SPM)
    observatory = EarthLocation(
            lat=31.0453*u.deg, lon=-115.4667*u.deg, height=2790*u.m)
    
    for t in time_1:
        
        # Alt/az reference frame at observatory
        frame_1 = AltAz(obstime=t, location=observatory)
    
        # Transform grid to alt/az coordinates at observatory
        altaz = radecs.transform_to(frame_1)
    
        strg[:,k] = (altaz.secz >= 0.) & (altaz.secz <= 3.)
        k += 1
    
    #print min(sum(strg,1)) 
    trg_slc = (sum(strg,1) < min_obsv/10.)
    #print sum(~trg_slc)
    ###
    if FIG:
        hpx, header = hp.read_map(File_Name, h=True, verbose=False)
	#
	plt.figure(1)
        hp.mollview(hpx)
        hp.projscatter(gwgc_ra_1,gwgc_dec_1,color='yellow',alpha=0.3,lonlat=True)
	hp.projscatter(RA[~trg_slc],DEC[~trg_slc],color='red',alpha=0.3,lonlat=True)
        plt.savefig(File_Name+'_'+str(nm)+'sigma.png')
    ###
    if (CLU):
    	fits.writeto('TEMP_'+File_Name+'_EM_Followup_'+str(nm)+'sigma'+'.fits',gwgc_slctd[1].data[gwgc_dec_allowed][wh_CLU][~trg_slc])
    else:
	fits.writeto('TEMP_'+File_Name+'_EM_Followup_'+str(nm)+'sigma'+'.fits',gwgc_slctd[1].data[gwgc_dec_allowed][~trg_slc])
    #    
    return tm_available, sum(~trg_slc) 		 
