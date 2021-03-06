import os, glob
import shutil
import healpy as hp
from numpy import sqrt,sum,where,pi,deg2rad,in1d,floor
from astropy.io import fits
from RATIR_EM_Followup_GW import GW_Followup 
from LIST_SLC_V2 import list_slc,list_slc_2  # LIST_SLC ---> LIST_SLC_V2
GO = True
GO2= False
CLU = True
def SKY_MAP_REDUCER(SKY_MAP,WHEN='',min_obsv=30.,obs_tm=5.,FIGs=False,Dist=1000.,REMOVE_TEMP=True):
	# comments:
	# min_obsv: should be a factor of 10 (10X)
	# obs_tm: observing time specified for each particular galaxy in the list [in minutes].
	# Dist: candidate galaxies within aLIGO region of Dist (Mpc) 
	#
	# the code uses current utc time unless you assign 
	# a specific time in WHEN='year-month-day hout-minute-second' format.
	# example: WHEN='2015-10-21 20:21:22'
	#
	hpx,header = hp.read_map(SKY_MAP, h=True, verbose=False)
	
	npix = len(hpx)
        print npix
	sky_area = 4 * 180**2 / pi
	#print "the sky resolution in arcmin: ", sqrt(sky_area / npix)*60.
	nside = hp.npix2nside(npix)

	################
	# GWGC Handeling
	################
	gwgc = fits.open('gwgc2011_2.fit')
	gwgc_ra = gwgc[1].data['RAJ2000']
	gwgc_dec = gwgc[1].data['DEJ2000']

	theta_c = 0.5 * pi - deg2rad(gwgc_dec)
	phi_c = deg2rad(gwgc_ra)
	ipix = hp.ang2pix(nside, theta_c, phi_c)
	################
	#
	j = hpx.argsort()[::-1]
	hpc = hpx[j].cumsum()
	fn = 1
	
	for sigma in [0.68,0.90,0.997]: # 0.95 ---> 0.90
    		# 
    		k = where(hpc<sigma)[0]
    		#print k.max()*1./len(hpx)
    		indx = j[:k.max()]
    		#	
   		ijp = in1d(ipix,indx)
		print sigma, sum(ijp),sky_area/npix * sum(ijp)
   		fits.writeto('TEMP_'+SKY_MAP+str(fn)+'sigma'+'.fits',gwgc[1].data[ijp])
		fn += 1
	#
	tm_avail, sz_lst = GW_Followup(SKY_MAP,WHEN=WHEN,min_obsv=min_obsv,nm=1,FIG=FIGs,CLU=CLU) # 1-sigma
	if ((floor(tm_avail*1./obs_tm) < sz_lst) and  not(GO)):
		shutil.copyfile('TEMP_'+SKY_MAP+'_EM_Followup_1sigma.fits','TEMP_'+SKY_MAP+'_EM_Followup_2sigma.fits')
		shutil.copyfile('TEMP_'+SKY_MAP+'_EM_Followup_1sigma.fits','TEMP_'+SKY_MAP+'_EM_Followup_3sigma.fits')
	else:
		tm_avail, sz_lst = GW_Followup(SKY_MAP,WHEN=WHEN,min_obsv=min_obsv,nm=2,FIG=FIGs,CLU=CLU) # 2-sigma
		if ((floor(tm_avail*1./obs_tm) < sz_lst) and not(GO)):
                	shutil.copyfile('TEMP_'+SKY_MAP+'_EM_Followup_2sigma.fits','TEMP_'+SKY_MAP+'_EM_Followup_3sigma.fits')
		else:
			tm_avail, sz_lst = GW_Followup(SKY_MAP,WHEN=WHEN,min_obsv=min_obsv,nm=3,FIG=FIGs,CLU=CLU) # 3-sigma
	#
	# tm_available: total available time for observations in minutes!
	# DARK_WINDOW would generate this and will pass it to the GW_Followup function!
	#print "\n total available observing time [min]: ", tm_avail
	tm_avail = 1.e6 # disable time_available concern
	if not(GO):
		list_slc(SKY_MAP,obs_tm=obs_tm,tm_available=tm_avail,Dist=Dist)
	else:
		list_slc_2(SKY_MAP,obs_tm=obs_tm,tm_available=tm_avail,Dist=Dist)
	#
	if REMOVE_TEMP:
		# removing all the Temporary files generated by the code!
		filelist = glob.glob("TEMP_*")
	        for f in filelist:
			os.remove(f)
