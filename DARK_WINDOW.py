# author: V. Zach Golkhou
from numpy import linspace,where,size,int,ceil,floor
from astropy import units as u
from astropy.time import Time
from datetime import datetime,timedelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun

def dark_window(tmrw=0,WHEN=''):
    # comments:
    # 
    #-------------------------------------------------
    if WHEN == '':
	t = datetime.utcnow()
    else:
	t = datetime.strptime(WHEN, "%Y-%m-%d %H:%M:%S")	
    time = Time(datetime(t.year, t.month, t.day+tmrw, 0, 0, 0)) #, scale='utc')

    delta_midnight = linspace(-400, 800, 121)*u.min
    time = time + delta_midnight
    
    # Geodetic coordinates of observatory (SPM)
    observatory = EarthLocation(
            lat=31.0453*u.deg, lon=-115.4667*u.deg, height=2790*u.m)
    
    # Alt/az reference frame at observatory
    frame = AltAz(obstime=time, location=observatory)
    
    # Where is the sun?
    sun_altaz = get_sun(time).transform_to(frame)
    
    wh = where((sun_altaz.alt <= -18*u.deg))[0]
    wh_sz = size(where(wh))
 
    night_strt = time[wh[0]]
    night_stop = time[wh[-1]]
    diff_1m = 0.0006944	
    tm_available = floor((night_stop.jd - night_strt.jd)*1./diff_1m) # in minutes	
    #print night_strt
    #print night_stop
    #
    if Time(t) > night_strt:
	#print "we've passed the night_start!"
	discard = timedelta(minutes=t.minute % 10,
                                 seconds=t.second,
                                 microseconds=t.microsecond)
    	t -= discard
    	if discard >= timedelta(minutes=5):
		t += timedelta(minutes=10)
       
    	t = Time(t)
	diff_10m = (t.jd - night_strt.jd)/0.007
	diff_10m = int(ceil(diff_10m)) 
	wh_sz = wh_sz - diff_10m
	night_strt = t
	if Time(t) < night_stop:
		tm_available = floor((night_stop.jd - t.jd)*1./diff_1m) # in minutes 

    return night_strt, wh_sz, tm_available    	
#-------------------------------------------------	
