import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
import scipy.optimize as sp
import math

def degreetorad(deg):
    '''
    This Procedure will convert Degree into Radian
    '''
    return deg*(np.pi / 180)


def mars_projections_on_ecliptic_plane(triangulation_data_index_pair,earth_holicentric_longitude,mars_holicentric_longitude):
    radiuslist=[]
    projectiondegreelist=[]
    position=[]
    for i in range(int(len(triangulation_data_index_pair)/2)):

        bita1 = earth_holicentric_longitude[i * 2]
        bita2 = earth_holicentric_longitude[i * 2 + 1]
        alpha1 = mars_holicentric_longitude[i * 2]
        alpha2 = mars_holicentric_longitude[i * 2 + 1]
        b1 = np.sin(degreetorad(bita1)) - np.tan(degreetorad(alpha1)) * np.cos(degreetorad(bita1))
        b2 = np.sin(degreetorad(bita2)) - np.tan(degreetorad(alpha2)) * np.cos(degreetorad(bita2))
        a2 = -1 * np.tan(degreetorad(alpha1))
        a4 = -1 * np.tan(degreetorad(alpha2))
        A = np.ones((2, 2))
        b = np.zeros((2, 1))
        b[0, 0] = b1
        b[1, 0] = b2
        A[0, 1] = a2
        A[1, 1] = a4
        coordinatex = np.matmul(np.linalg.inv(A), b)[1, 0]
        coordinatey = np.matmul(np.linalg.inv(A), b)[0, 0]
        position.append([coordinatex,coordinatey])
        r = np.sqrt(np.square(coordinatex) + np.square(coordinatey))
        if coordinatex > 0 and coordinatey > 0:
            phi = np.arctan(coordinatey / coordinatex)
        elif coordinatex > 0 > coordinatey:
            phi = 2 * np.pi - np.arctan(np.abs(coordinatey) / coordinatex)
        elif coordinatey > 0 > coordinatex:
            phi = np.pi - np.arctan(coordinatey / np.abs(coordinatex))
        else:
            phi = np.pi - np.arctan(np.abs(coordinatey) / np.abs(coordinatex))
        projectiondegree = (180 / np.pi) * phi
        radiuslist.append(r)
        projectiondegreelist.append(projectiondegree)
    return radiuslist,projectiondegreelist,position

def fitcircle(radiuslist,radius):
    x0=[radius]
    params = sp.minimize(func_calculate_squared_loss,x0,args=[radiuslist], bounds=[(1, None)])
    return params['x'],params['fun']

def func_calculate_squared_loss(params,args):
    loss=0
    radius=params[0]
    radiuslist=args[0]
    for i in range(len(radiuslist)):
        loss=loss+math.pow(radius - radiuslist[i], 2)
    return loss


if __name__ == "__main__":
    mars_opposition_data = pd.read_csv('./../data/01_data_mars_opposition.csv')
    mars_triangulation_data = pd.read_csv('./../data/01_data_mars_triangulation.csv')

    triangulation_data_index_pair = mars_triangulation_data['PairIndex'].values

    earth_HelioCentric_degree = mars_triangulation_data['DegreeEarthLocationHelioCentric'].values
    earth_HelioCentric_min = mars_triangulation_data['MinuteEarthLocationHelioCentric'].values
    earth_holicentric_longitude=earth_HelioCentric_degree+(earth_HelioCentric_min/60)

    mars_HelioCentric_degree = mars_triangulation_data['DegreeMarsLocationGeoCentric'].values
    mars_HelioCentric_min = mars_triangulation_data['MinuteMarsLocationGeoCentric'].values
    mars_holicentric_longitude=mars_HelioCentric_degree+(mars_HelioCentric_min/60)

    radiuslist, projectiondegreelist, position=mars_projections_on_ecliptic_plane(triangulation_data_index_pair,
                                                                                  earth_holicentric_longitude,
                                                                                  mars_holicentric_longitude)
    print("2(i) Different projections of Mars's location on the ecliptic plane:[X,Y]")
    print("")
    print(position)
    print("")
    optimised_orbit_radius,loss=fitcircle(radiuslist,1)
    print("")
    print("2(ii) Optimised Radius(in AU) For best fit circle centred at the Sun.")
    print("")
    print(optimised_orbit_radius)
    print("")
    print("Loss function :Squared Euclidean distance  ")
    print("2(ii) Total Loss."+str(loss))