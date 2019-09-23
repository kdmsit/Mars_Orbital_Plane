import numpy as np
import pandas as pd
import scipy.optimize as sp
import matplotlib.patches as ptch
import matplotlib.pyplot as plt

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

def calculate_zcoord(xcoord_list,ycoord_list,a,b,c):
    zcoord_list=[]
    for i in range(len(xcoord_list)):
        zcoord=(-a*xcoord_list[i]-b*ycoord_list[i])/c
        zcoord_list.append(zcoord)
    return zcoord_list

def compute_loss(params,args):
    xlist, ylist=args[0],args[1]
    foci2x = np.asarray(params[0])
    foci2y = np.asarray(params[1])
    majoraxix=np.asarray(params[2])
    foci2 = np.asarray([foci2x, foci2y])
    loss=[]
    for i in range(5):
        a=np.sqrt(np.square(xlist[i])+np.square(ylist[i]))
        b=np.sqrt(np.square(xlist[i]-foci2x)+np.square(ylist[i]-foci2y))
        loss.append(abs((a+b)-majoraxix))
    totalloss= np.sum(loss)
    #print(totalloss)
    return totalloss
def min_fit_ellipse(xcoord_list,ycoord_list,foci2x,foci2y,majoraxix):

    x0 = [foci2x,foci2y,majoraxix]
    params = sp.minimize(compute_loss, x0,method="L-BFGS-B",args=[xcoord_list, ycoord_list])
    return params['x'],params['fun']

def find_best_fit_ellipse_on_mars_orbital_plane(xcoord_list,ycoord_list):
    foci2x = 1
    foci2y = 1
    majoraxix=1
    return min_fit_ellipse(xcoord_list,ycoord_list,foci2x,foci2y,majoraxix)


if __name__ == "__main__":
    # region All Datas
    inclination_angle_deg=1.8536
    # Opposition Coordinates found in Question 1.
    opp_ylist = [1.4394458945438964,1.5014896753783606,0.9747386019545708,0.11722917382644234,-0.8865956950779831,
            -1.566917999583486,-0.47785030102239034,1.1580390999650196,1.5681475650401604,1.2259167541156522,
            0.4732490965051435,-0.5012114194578124]
    opp_xlist=[0.6265943464442608,-0.4569031621546369,-1.2299365359470285,-1.5651874087757587,-1.2956486880860487,
               -0.08989061772290922,1.4943008401709512,1.0601149388047009,-0.06755277403225604,-0.9798047440522887,
            -1.4963997424804596,-1.4876493315653687]
    # Triangulation Coordinates and radius found in Question 2.
    xcoord_list=[-1.4529736727603795,1.195672782788594,1.073885314206997,-1.6323045900130564,-1.5537673314861349]
    ycoord_list=[0.865533530153104,-0.6868566346181091,1.0511069275483509,-0.14854179871578332,0.6248989852957587]
    radius_list=[1.6912364664806872, 1.3789145876450561, 1.502682748024211, 1.6390494014957, 1.674721368531856]
    a=0.04709835186731311
    b=-0.04404339492276171
    c=1.9924440262398444
    radius = 1.57732091
    loss=0.07120993396639318
    # endregion

    # region Question 4(i)
    zcoord_list=calculate_zcoord(xcoord_list,ycoord_list,a,b,c)
    print("4(i) Mars's five different (3-d) locations on Mars's orbital plane :")
    print("X_Coord :" + str(xcoord_list))
    print("Y_Coord :" + str(ycoord_list))
    print("Z_Coord :" + str(zcoord_list))
    print("")
    # endregion

    # region Question 4(ii)
    print("4(ii) Best fit circle on Mars's orbital plane with the Sun as centre :")
    print("Radius :" + str(round(radius,2)))
    print("Total Loss :" + str(round(loss,4)))
    print("")
    # endregion

    # region Question 4(iii)
    opt_parameters,total_loss=find_best_fit_ellipse_on_mars_orbital_plane(xcoord_list,ycoord_list)
    foci2x = opt_parameters[0]
    foci2y = opt_parameters[1]
    print("4_(iii) best fit ellipse on Mars's orbital plane is :")
    print(" Sum of total loss :")
    print(total_loss)
    majoraxix = opt_parameters[2]
    c=np.sqrt(np.square(foci2x)+np.square(foci2y))/2
    height=2*np.sqrt(np.square(majoraxix/2)-np.square(c/2))
    print("Foci1 as Sun: " + str([0,0]) + " ,Other Focai: " + str([foci2x,foci2y]))
    print("Major Axis: "+str(majoraxix)+" Minor Axis: "+str(height))
    print("")
    # endregion

    # region Plots
    ellipse=ptch.Ellipse((foci2x/2,foci2y/2), majoraxix, height, angle=0,color="r",fill=False,label='Ellipse')
    circle1 = plt.Circle((0, 0), radius, color='b', fill=False,label='Circle')
    ax = plt.gca()
    ax.cla()
    ax.add_artist(ellipse)
    ax.add_artist(circle1)
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    plt.grid(linestyle='--')
    for i in range(len(xcoord_list)):
        if(i==0):
            plt.scatter(xcoord_list[i], ycoord_list[i], s=20, facecolors='b', edgecolors='b',label='triangulation mars position')
        else:
            plt.scatter(xcoord_list[i], ycoord_list[i], s=20, facecolors='b', edgecolors='b')
    plt.legend()
    for i in range(len(opp_xlist)):
        if(i==0):
            plt.scatter(opp_xlist[i], opp_ylist[i], s=20, facecolors='y', edgecolors='g',label='Opposition mars position')
        else:
            plt.scatter(opp_xlist[i], opp_ylist[i], s=20, facecolors='y', edgecolors='y')
    plt.legend()
    plt.scatter(0, 0, s=50, facecolors='r', edgecolors='r',label='Sun position')
    plt.legend()
    plt.title('Plot of 4(ii) and 4(iii) with the Sun at one of the foci[Dark Red at (0,0)].', fontsize=8)
    plt.savefig("4.png", bbox_inches='tight')
    plt.show()
    print("Plot of 4(ii) and 4(iii) saved in home directory")
    print("")
    plt.close()
    # endregion
