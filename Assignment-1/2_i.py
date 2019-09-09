import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
import scipy.optimize as sp

def degreetorad(deg):
    '''
    This Procedure will convert Degree into Radian
    '''
    return deg* (np.pi / 180)

def calculate_longitude(zodiac_indexes,degree,minute,second):
    '''
    This  Procedure will calculate Longitude given zodiac_indexes,degree,minute,second
    '''
    return zodiac_indexes * 30 + degree + (
            minute / 60) + (second / 3600)

def diff_logarithmean_loggeomean_radius(parameters,args):
    '''
    This Procedure will calculate radius for each point of mars and calculate Loss[log(arithmetic_mean) - log(geometric_mean)]
    '''
    x = parameters[0]
    alpha = args[0]
    bita = args[1]
    y = parameters[1]
    radius_list = []
    for i in range(len(alpha)):
        z1 = x * np.sin(degreetorad(bita[i] - y))
        z2 = np.sin(degreetorad(alpha[i] - y))
        z3 = np.cos(degreetorad(bita[i] - alpha[i]))
        z4 = 2 * z1 * z2 * z3 + np.square(z1) + np.square(z2)
        z5 = 1 - np.square(z3)
        radius = np.sqrt(z4 / z5)
        radius_list.append(radius)
    arithmetic_mean = np.mean(radius_list)
    geometric_mean = gmean(radius_list)
    return np.log(arithmetic_mean) - np.log(geometric_mean)

def find_radius(alpha,bita,x_opt,y_opt):
    '''
    Given Optimised x and y this procedure will find the optimum radius value.
    '''
    radius_list=[]
    for i in range(12):
        rsin1 = x_opt * np.sin(degreetorad(bita[i] - y_opt))
        rsin2 = np.sin(degreetorad(alpha[i] - y_opt))
        cos3 = np.cos(degreetorad(bita[i] - alpha[i]))
        z1 = 2 * rsin1 * rsin2 * cos3 + np.square(rsin1) + np.square(rsin2)
        z2 = 1 - np.square(cos3)
        radius = np.sqrt(z1 / z2)
        radius_list.append(radius)
    print(radius_list)

def minimize_loss(alpha,bita,x,y)  :
    '''
    Given Optimised x and y this procedure will find the optimum radius value.
    This is just to confirm that all radiuses is of approximately same length.
    '''
    parameters=[x,y]
    optimised_res = sp.minimize(diff_logarithmean_loggeomean_radius, parameters,args=[alpha, bita],bounds=((0,None), (None,None)))
    return optimised_res.x[0],optimised_res.x[1],optimised_res.fun

if __name__ == "__main__":
    mars_data = pd.read_csv('./../data/01_data_mars_opposition.csv')

    # region Calculate Longitude Relative to Actual Sun(ALPHA)
    actual_sun_zodiac_indexes=mars_data['ZodiacIndex'].values
    actual_sun_degree = mars_data['Degree'].values
    actual_sun_minute = mars_data['Minute'].values
    actual_sun_second = mars_data['Second'].values
    longitude_relative_actual_sun=calculate_longitude(actual_sun_zodiac_indexes,actual_sun_degree,actual_sun_minute,actual_sun_second)
    # endregion

    # region Calculate Longitude Relative to Average Sun(BITA)
    avg_sun_zodiac_indexes = mars_data['ZodiacIndexAverageSun'].values
    avg_sun_degree = mars_data['DegreeMean'].values
    avg_sun_minute = mars_data['MinuteMean'].values
    avg_sun_second = mars_data['SecondMean'].values
    longitude_relative_avg_sun = calculate_longitude(avg_sun_zodiac_indexes, avg_sun_degree, avg_sun_minute,avg_sun_second)
    # endregion
    #Initial Guess of x :1.2, y :120
    x=1.2
    y=120
    print("Initial Guess of x :"+str(x)+", y :"+str(y))
    x_opt, y_opt,loss=minimize_loss(longitude_relative_actual_sun,longitude_relative_avg_sun,x,y)
    print("Optimised Value of x :" + str(x_opt) + ", y :" + str(y_opt))
    print("Total Loss :"+str(loss))
    #ANSWRE: Optimised Value of x :0.9662028648251194, y :148.874455333893
    #find_radius(longitude_relative_actual_sun, longitude_relative_avg_sun, x_opt, y_opt)