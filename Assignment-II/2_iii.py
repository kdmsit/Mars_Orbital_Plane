import numpy as np
import pandas as pd
import scipy.optimize as sp

def degreetorad(deg):
    '''
    This Procedure will convert Degree into Radian
    '''
    return deg* (np.pi / 180)

def radtodegree(rad):
    '''
    This Procedure will convert Radian into Degree
    '''
    return rad * (180 / np.pi)


def calculate_longitude(zodiac_indexes,degree,minute,second):
    '''
    This  Procedure will calculate Longitude given zodiac_indexes,degree,minute,second
    '''
    return zodiac_indexes * 30 + degree + (
            minute / 60) + (second / 3600)

def find_mars_holicentric_latitude(mars_geocentric_latitude,radius):
    mars_holicentric_latitude=[]
    for latitude in mars_geocentric_latitude:
        latitude=degreetorad(latitude)
        tanthita=np.tan(latitude)
        r=(radius-1)/radius
        holicentric_latitude=np.arctan(r*tanthita)
        mars_holicentric_latitude.append(radtodegree(holicentric_latitude))
    return mars_holicentric_latitude

def find_projetion(mars_holicentric_latitude,longitude_relative_actual_sun):
    xlist = []
    ylist = []
    zlist = []
    for i in range(len(longitude_relative_actual_sun)):
        x = radius * np.cos(degreetorad(mars_holicentric_latitude[i])) * np.cos(
            degreetorad(longitude_relative_actual_sun[i]))
        xlist.append(x)
        y = radius * np.cos(degreetorad(mars_holicentric_latitude[i])) * np.sin(
            degreetorad(longitude_relative_actual_sun[i]))
        ylist.append(y)
        z = radius * np.sin(degreetorad(mars_holicentric_latitude[i]))
        zlist.append(z)
    return xlist,ylist,zlist

def compute_total_distance(params,args):        #xlist,ylist,zlist
    a = params[0]
    b = params[1]
    c = params[2]
    xlist, ylist, zlist=args[0],args[1],args[2]
    distance=[]
    for i in range(len(xlist)):
        nu=a*xlist[i]+b*ylist[i]+c*zlist[i]
        de=np.sqrt(a**2+b**2+c**2)
        distance.append(np.square(abs(nu)/de))
    return np.sum(distance)

def compute_optimised_parameters(xlist, ylist, zlist):
    a = 1
    b = 1
    c = 1
    x0 = [a,b,c]
    params = sp.minimize(compute_total_distance, x0, args=[xlist,ylist,zlist])
    return params['x']

if __name__ == "__main__":
    mars_opposition_data = pd.read_csv('./../data/01_data_mars_opposition.csv')
    mars_triangulation_data = pd.read_csv('./../data/01_data_mars_triangulation.csv')

    # region Read from Mars Opposition Data
    actual_sun_zodiac_indexes = mars_opposition_data['ZodiacIndex'].values
    actual_sun_degree = mars_opposition_data['Degree'].values
    actual_sun_minute = mars_opposition_data['Minute'].values
    actual_sun_second = mars_opposition_data['Second'].values
    longitude_relative_actual_sun = calculate_longitude(actual_sun_zodiac_indexes, actual_sun_degree, actual_sun_minute,
                                                        actual_sun_second)
    actual_sun_lat_deg=mars_opposition_data['LatDegree'].values
    actual_sun_lat_min = mars_opposition_data['LatMinute'].values
    mars_geocentric_latitude = actual_sun_lat_deg + (actual_sun_lat_min / 60)
    # endregion

    # region Read From Mars Triangulation Data
    triangulation_data_index_pair = mars_triangulation_data['PairIndex'].values
    earth_HelioCentric_degree = mars_triangulation_data['DegreeEarthLocationHelioCentric'].values
    earth_HelioCentric_min = mars_triangulation_data['MinuteEarthLocationHelioCentric'].values
    earth_holicentric_longitude=earth_HelioCentric_degree+(earth_HelioCentric_min/60)
    mars_HelioCentric_degree = mars_triangulation_data['DegreeMarsLocationGeoCentric'].values
    mars_HelioCentric_min = mars_triangulation_data['MinuteMarsLocationGeoCentric'].values
    mars_holicentric_longitude=mars_HelioCentric_degree+(mars_HelioCentric_min/60)
    # endregion

    radius = 1.57732091

    # region 3(i) Heliocentric Latitudes of Mars
    mars_holicentric_latitude=find_mars_holicentric_latitude(mars_geocentric_latitude,radius)
    print("3(i) 12 Heliocentric Latitudes of Mars.")
    print("")
    for latitude in mars_holicentric_latitude:
        print(latitude)
    print("")
    # endregion
    xlist, ylist, zlist=find_projetion(mars_holicentric_latitude,longitude_relative_actual_sun)
    opt_parameters=compute_optimised_parameters(xlist, ylist, zlist)
    print("3(ii) Best fit for Mars's orbital plane is ax+by+cz=0 where optimised parameters a,b,c are :.")
    print("a : "+str(opt_parameters[0]))
    print("b : " + str(opt_parameters[1]))
    print("c : " + str(opt_parameters[2]))
    print("")
    inclination=np.arccos(opt_parameters[2]/np.linalg.norm(opt_parameters))
    inclination=radtodegree(inclination)
    inclinationdegree=int(inclination)
    inclinationmin=(inclination - inclinationdegree)*60
    print("3(iii) Inclination of this plane in degrees/minutes:")
    print(str(inclinationdegree)+" degree")
    print(str(inclinationmin) + " minutes")