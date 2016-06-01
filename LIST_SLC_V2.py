from pyfits import getdata
from numpy import array,arange,size,savetxt,where,max,sum,argsort,floor

dist_sort = False
bmag_sort = True

def list_slc(FILE_NAME,obs_tm=5,tm_available=10.,Dist=10.):
	"""
	obs_tm: observing time for each specific galaxy in minutes.
	tm_available: total available time for observations in minutes!
	Dist: candidate galaxies within aLIGO region of Dist (Mpc) 
	"""
	data_1,hdr_1 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'1sigma'+'.fits',header=True)
        data_2,hdr_2 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'2sigma'+'.fits',header=True)
        data_3,hdr_3 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'3sigma'+'.fits',header=True)
        #
	Dist_1 = data_1['Dist']
	Dist_2 = data_2['Dist']
	Dist_3 = data_3['Dist']	
	#
	wh_1 = (Dist_1 <= Dist)
	wh_2 = (Dist_2 <= Dist)
	wh_3 = (Dist_3 <= Dist)
	#
	Ngal = floor(tm_available*1./obs_tm) # the number of galaxies that we can observe for the time_available
	#print sum(wh_1), sum(wh_2), sum(wh_3)
	#print "Time Available (min): ", tm_available
	#sz_1 = size(data_1[wh_1])
	#sz_2 = size(data_2[wh_2])
	#sz_3 = size(data_3[wh_3])
	#
	# sz_1 <= sz_2 <= sz_3
	gs = array([sum(wh_1),sum(wh_2),sum(wh_3)])*1.*obs_tm/tm_available
	if ((gs[0] > 1) | (gs[0]==gs[1]==gs[2])):
		pck = 1
	else:
		pck = max(where(gs<=1)[0])+1
		if (gs[1]==gs[2]):
			pck = 2
	
        if pck == 1: 
		RA_slc = data_1['RAJ2000'][wh_1]; DEC_slc = data_1['DEJ2000'][wh_1]; Dist_slc = Dist_1[wh_1]; BMAG_slc = data_1['BMAG1'][wh_1]
	elif pck == 2:
		RA_slc = data_2['RAJ2000'][wh_2]; DEC_slc = data_2['DEJ2000'][wh_2]; Dist_slc = Dist_2[wh_2]; BMAG_slc = data_2['BMAG1'][wh_2]
        else:
		RA_slc = data_3['RAJ2000'][wh_3]; DEC_slc = data_3['DEJ2000'][wh_3]; Dist_slc = Dist_3[wh_3]; BMAG_slc = data_3['BMAG1'][wh_3]
        #
        #wh = (Dist_slc <= Dist)
	if dist_sort: ags = argsort(Dist_slc)
	if bmag_sort: ags = argsort(BMAG_slc)
        if (Ngal < size(RA_slc)):
		chs = Ngal
	else:
		chs = size(RA_slc)
	# 
	RA_DEC_slc = [[RA_slc[ags][:chs][k],DEC_slc[ags][:chs][k],Dist_slc[ags][:chs][k],BMAG_slc[ags][:chs][k]] for k in arange(size(RA_slc[:chs]))]
        #
        savetxt('galx_trg_lst_'+str(pck)+'sigma_'+str(int(chs))+'Ngalx.txt',RA_DEC_slc,fmt='%.4f',newline='\n')

def list_slc_2(FILE_NAME,obs_tm=5,tm_available=10.,Dist=10.):
        """
        obs_tm: observing time for each specific galaxy in minutes.
        tm_available: total available time for observations in minutes!
        Dist: candidate galaxies within aLIGO region of Dist (Mpc) 
        """
        data_1,hdr_1 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'1sigma'+'.fits',header=True)
        data_2,hdr_2 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'2sigma'+'.fits',header=True)
        data_3,hdr_3 = getdata('TEMP_'+FILE_NAME+'_EM_Followup_'+'3sigma'+'.fits',header=True)
        #
        Dist_1 = data_1['Dist']
        Dist_2 = data_2['Dist']
        Dist_3 = data_3['Dist']
        #
        wh_1 = (Dist_1 <= Dist)
        wh_2 = (Dist_2 <= Dist)
        wh_3 = (Dist_3 <= Dist)
        #
        Ngal = floor(tm_available*1./obs_tm) # the number of galaxies that we can observe for the time_available
        
        RA_slc_1 = data_1['RAJ2000'][wh_1]; DEC_slc_1 = data_1['DEJ2000'][wh_1]; Dist_slc_1 = Dist_1[wh_1]; BMAG_slc_1 = data_1['BMAG1'][wh_1]
        RA_slc_2 = data_2['RAJ2000'][wh_2]; DEC_slc_2 = data_2['DEJ2000'][wh_2]; Dist_slc_2 = Dist_2[wh_2]; BMAG_slc_2 = data_2['BMAG1'][wh_2]
        RA_slc_3 = data_3['RAJ2000'][wh_3]; DEC_slc_3 = data_3['DEJ2000'][wh_3]; Dist_slc_3 = Dist_3[wh_3]; BMAG_slc_3 = data_3['BMAG1'][wh_3]
        #
        #wh = (Dist_slc <= Dist)
	if dist_sort: ags_1 = argsort(Dist_slc_1)
	if bmag_sort: ags_1 = argsort(BMAG_slc_1)
        if (Ngal < size(RA_slc_1)):
                chs = Ngal
        else:
                chs = size(RA_slc_1)

        RA_DEC_slc_1 = [[RA_slc_1[ags_1][:chs][k],DEC_slc_1[ags_1][:chs][k],Dist_slc_1[ags_1][:chs][k],BMAG_slc_1[ags_1][:chs][k]] for k in arange(size(RA_slc_1[:chs]))]
        savetxt('galx_trg_lst_'+str(1)+'sigma_'+str(int(chs))+'Ngalx.txt',RA_DEC_slc_1,fmt='%.4f',newline='\n')
	#
        if dist_sort: ags_2 = argsort(Dist_slc_2)
	if bmag_sort: ags_2 = argsort(BMAG_slc_2)
        if (Ngal < size(RA_slc_2)):
                chs = Ngal
        else:
                chs = size(RA_slc_2)
        
        RA_DEC_slc_2 = [[RA_slc_2[ags_2][:chs][k],DEC_slc_2[ags_2][:chs][k],Dist_slc_2[ags_2][:chs][k],BMAG_slc_2[ags_2][:chs][k]] for k in arange(size(RA_slc_2[:chs]))]
        savetxt('galx_trg_lst_'+str(2)+'sigma_'+str(int(chs))+'Ngalx.txt',RA_DEC_slc_2,fmt='%.4f',newline='\n')	
	#
        if dist_sort: ags_3 = argsort(Dist_slc_3)
	if bmag_sort: ags_3 = argsort(BMAG_slc_3)
        if (Ngal < size(RA_slc_3)):
                chs = Ngal
        else:
                chs = size(RA_slc_3)

        RA_DEC_slc_3 = [[RA_slc_3[ags_3][:chs][k],DEC_slc_3[ags_3][:chs][k],Dist_slc_3[ags_3][:chs][k],BMAG_slc_3[ags_3][:chs][k]] for k in arange(size(RA_slc_3[:chs]))]
        savetxt('galx_trg_lst_'+str(3)+'sigma_'+str(int(chs))+'Ngalx.txt',RA_DEC_slc_3,fmt='%.4f',newline='\n')
        #
