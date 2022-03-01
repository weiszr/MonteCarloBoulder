from __future__ import print_function
from numpy import *
from random import uniform
import tqdm
from random import  uniform
#import sys
import time
from multiprocessing import Pool
from numpy import linspace,sqrt,zeros
from tqdm import tqdm_notebook
from random import  uniform
from multiprocessing import Pool
import matplotlib.pyplot as plt
from src.helpers import *

    
#--------------------------------------Jen--------------------------------------------
##Rotate coordinates
def rotatecoord(xin,yin,th): 
#	from numpy import sin, cos
    xout = xin*cos(th) + yin*sin(th)
    yout = -xin*sin(th) + yin*cos(th)
    return xout,yout


def trapezoidchar(aa,bb,cc,hh):
    """
    Calculate area and centroid of trapezoid
    """
    A = hh*(aa + bb)/2.
    xc = (2.*aa*cc + aa**2. + cc*bb + aa*bb + bb**2.)/3./(aa + bb)
    yc = hh*(2.*aa + bb)/3./(aa + bb)
    return A,xc,yc

def blockcorners(b,c,lo,th,thetacr):
	"""
	Calculate corners of block
	"""
	from numpy import sin, cos
#	from numpy import sin, cos
	x2 = -b*cos(th)
	z2 = b*sin(th)
	x3 = -2.*lo*sin(thetacr - th)
	z3 = 2.*lo*cos(thetacr - th)
	x4 = c*sin(th)
	z4 = c*cos(th)
	return x2,z2,x3,z3,x4,z4


def shapemoments(b,c,lo,eta,th,thetacr):
    """
    Calculate area, centrid, and moment of inertia for submerged portion of block
    using python polygon functions. Function calculates quantities only when th is
    at or just over thetacr, zeros are returned for other th. Function is coded for 
    quadrilaterals whose lower-right corner is at (x,z) = (0,0). The polygon portion
    of the code is modified from:
    armatita at: https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
    Error statistics for a 6x6 m block cross-section based on nresx
    were estimated as follows:
        nresx=5000., A: 0.03% error, areaI: 0.06% error
        nresx=1000., A: 0.06% error, areaI: 0.12% error
        nresx=500.,  A: 0.12% error, areaI: 0.22% error, 0.22% error in wettedx    
        nresx=250.,  A: 1.03% error, areaI: 1.84% error, 0.40% error in wettedx       
    """
    from numpy import linspace,meshgrid,newaxis,ones,hstack,shape,sum,add,power,max,min,pi,tan
    from matplotlib import path
    if th <= (thetacr+pi/100.):
        nresx = 500. 
        nresz = int(lo*nresx/(b+c))
        x = linspace(-b,c,nresx)
        z = linspace(0.,2*lo,nresz)
        dx = x[1] - x[0]
        dz = z[1] - z[0]
        x = x + 0.5*dx
        z = z + 0.5*dx
        xv0,zv0 = meshgrid(x,z)
        xv0 = xv0.flatten()
        zv0 = zv0.flatten()
        xv = xv0[zv0<=eta][:,newaxis]
        zv = zv0[zv0<=eta][:,newaxis]
        dA = dx*dz*ones(shape(xv)[0])
        xdA = xv*dx*dz
        zdA = zv*dx*dz
        r2dA = dx*dz*add(power(xv,2.),power(zv,2))
        x2,z2,x3,z3,x4,z4 = blockcorners(b,c,lo,th,thetacr)
        p = path.Path([(0.,0.), (x2,z2), (x3,z3), (x4,z4)]) 
        flags = p.contains_points(hstack((xv,zv)))
        A = sum(dA[flags])
        xc = sum(xdA[flags])/A
        zc = sum(zdA[flags])/A
        areaI = sum(r2dA[flags]) 
        wettedx = max(xv[flags]) - min(xv[flags]) + dx
        xtil = (max(xv[flags]) + min(xv[flags]))/2.
    else:
        A = 0.
        xc = 0.
        zc = 0.
        areaI = 0.
        wettedx = 0.
        xtil = 0.
    return A,xc,zc,areaI,wettedx,xtil


#Calculation of area, centroid, and moment of inertia using block side slopes
def blockmoments(b,c,lo,eta,th,thetacr):
    """
    Calculate area, centroid, and moment of inertia for submerged portion of block using 
    cross-section sides. Function calculates quantities only when th is at or just over 
    thetacr, zeros are returned for other th. Function is coded for rectangles whose
    lower-right corner is at (x,z) = (0,0). Error statistics for a 6x6 m block 
    cross-section based on nresx were estimated as follows:
        nresx=5000., A: 0.03% error, areaI: 0.06% error
        nresx=1000., A: 0.06% error, areaI: 0.12% error
        nresx=500.,  A: 0.12% error, areaI: 0.22% error, 0.22% error in wettedx    
        nresx=250.,  A: 1.03% error, areaI: 1.84% error, 0.40% error in wettedx       
    """
    from numpy import linspace,meshgrid,ones,shape,sum,add,power,max,min,pi,tan
    if th <= (thetacr+pi/100.): 
        nresx = 500. 
        nresz = int(lo*nresx/(b+c))
        x = linspace(-b,c,nresx)
        z = linspace(0.,2*lo,nresz)
        dx = x[1] - x[0]
        dz = z[1] - z[0]
        x = x + 0.5*dx
        z = z + 0.5*dx
        xv0,zv0 = meshgrid(x,z)
        xv0 = xv0.flatten()
        zv0 = zv0.flatten()
        xv = xv0[zv0<=eta]
        zv = zv0[zv0<=eta]
        x2,z2,x3,z3,x4,z4 = blockcorners(b,c,lo,th,thetacr)
        m23 = tan(pi/2.-th)
        b23 = z2 - m23*x2
        xv1 = xv[zv-(m23*xv+b23) <= 0.]
        zv1 = zv[zv-(m23*xv+b23) <= 0.]
        m14 = tan(pi/2.-th)
        b14 = 0.
        xv2 = xv1[zv1-(m14*xv1+b14) >= 0.]
        zv2 = zv1[zv1-(m14*xv1+b14) >= 0.]
        m12 = tan(pi-th)
        b12 = 0.
        xv3 = xv2[zv2-(m12*xv2+b12) >= 0.]
        zv3 = zv2[zv2-(m12*xv2+b12) >= 0.]
        m43 = tan(pi-th)
        b43 = z4 - m43*x4
        xv4 = xv3[zv3-(m43*xv3+b43) <= 0.]
        zv4 = zv3[zv3-(m43*xv3+b43) <= 0.]
#debug         if (th > 0.788) and (th < 0.789): debugplots(x2,x3,x4,z2,z3,z4,xv,zv,xv1,zv1,xv2,zv2,xv3,zv3,xv4,zv4)
        dA = dx*dz*ones(shape(xv4)[0])
        xdA = xv4*dx*dz
        zdA = zv4*dx*dz
        r2dA = dx*dz*add(power(xv4,2.),power(zv4,2))
        A = sum(dA)
        xc = sum(xdA)/A
        zc = sum(zdA)/A
        areaI = sum(r2dA)
        wettedx = max(xv4) - min(xv4) + dx
        xtil = (max(xv4) + min(xv4))/2.
    else:
        A = 0.
        xc = 0.
        zc = 0.
        areaI = 0.
        wettedx = 0.
        xtil = 0.
    return A,xc,zc,areaI,wettedx,xtil

#JInote: need to pass in rotated coordinates with origin at 0 and -1 indices
#should check th<=thetacr before calling function
def polymoments(xrot,zrot,lo,eta):
    """
    Calculate area, centroid, and (mass moment of inertia/prism width/density) for 
    submerged portion of polygon prism by breaking area into right triangle areas.
    Function assumes the vector angle from (0,0) to two adjacent 
    verticies is <= 90 deg. The 1D arrays xrot and zrot are the polygon vertices, where
    the first vertice in lower right position and appearing twice, at index=0 and index=-1. 
    Error statistics are:
    """
    from numpy import array,cos,sin,tan,pi,zeros,argsort,argwhere,insert,delete,min,max
    from math import acos,asin
    iabove = argwhere(zrot>eta)
    if len(iabove)>0:
        zeta1 = eta
        m1 = (xrot[min(iabove)]-xrot[min(iabove)-1])/(zrot[min(iabove)]-zrot[min(iabove)-1])
        xeta1 = (eta-zrot[min(iabove)-1])*m1 + xrot[min(iabove)-1]
        zeta2 = eta 
        m2 = (xrot[max(iabove)+1]-xrot[max(iabove)])/(zrot[max(iabove)+1]-zrot[max(iabove)])
        xeta2 = (eta-zrot[max(iabove)])*m2 + xrot[max(iabove)]
        x = delete(xrot,iabove)
        z = delete(zrot,iabove)
        x = insert(x,min(iabove),[xeta1,xeta2])
        z = insert(z,min(iabove),[zeta1,zeta2])
    else:
        x = xrot
        z = zrot
    MOI = 0.
    TotalArea = 0.
    cxSum = -sum( (x[0:-1]+x[1::])*(x[0:-1]*z[1::]-x[1::]*z[0:-1]) ) #negative due to summation order
    czSum = -sum( (z[0:-1]+z[1::])*(x[0:-1]*z[1::]-x[1::]*z[0:-1]) ) #negative due to summation order
    for ii in range(len(x)-3):
        i = ii + 1
        leng = zeros(3)
        ang = zeros(3)
        xT = array([0.,x[i],x[i+1]])
        zT = array([0.,z[i],z[i+1]])
        leng[2] = (xT[1]**2. + zT[1]**2.)**0.5
        leng[1] = (xT[2]**2. + zT[2]**2.)**0.5
        leng[0] = ( (xT[2]-xT[1])**2. + (zT[2]-zT[1])**2. )**0.5
        dotprod = xT[1]*xT[2] + zT[1]*zT[2]
        arg0 = dotprod/leng[2]/leng[1]
        ang[0] = acos(max([min([arg0,1.]) , -1.]) ) #to avoid "math domain error"    
        arg1 = leng[1]*sin(ang[0])/leng[0]
        if leng[1]<=leng[2]:
            arg0 = leng[1]*sin(ang[0])/leng[0]
            ang[1] = asin(max([min([arg0,1.]) , -1.]) ) 
            ang[2] = pi - ang[0] - ang[1] #asin(leng[2]*sin(ang[0])/leng[0])
        else:
            arg0 = leng[2]*sin(ang[0])/leng[0]
            ang[2] = asin(max([min([arg0,1.]) , -1.]) ) 
            ang[1] = pi - ang[0] - ang[2] 
        isort = argsort(leng)[::-1]
        xT = xT[isort]
        zT = zT[isort]
        leng = leng[isort]
        ang = ang[isort]
        L1 = leng[0]/(1./tan(ang[1]) + 1./tan(ang[2]))
        L2_T1 = L1/tan(ang[1])
        L2_T2 = L1/tan(ang[2])
        xo = (L2_T1/leng[0])*(xT[2]-xT[1]) + xT[1]
        zo = (L2_T1/leng[0])*(zT[2]-zT[1]) + zT[1]
        xc_T1 = (xo+xT[0]+xT[1])/3.
        zc_T1 = (zo+zT[0]+zT[1])/3.
        xc_T2 = (xo+xT[0]+xT[2])/3.
        zc_T2 = (zo+zT[0]+zT[2])/3.
        area_T1 = 0.5*L1*L2_T1
        area_T2 = 0.5*L1*L2_T2
        I_T1 = area_T1*(L1**2 + L2_T1**2)/6. - ((xo-xc_T1)**2.+(zo-zc_T1)**2.)*area_T1 + (xc_T1**2. + zc_T1**2.)*area_T1
        I_T2 = area_T2*(L1**2 + L2_T2**2)/6. - ((xo-xc_T2)**2.+(zo-zc_T2)**2.)*area_T2 + (xc_T2**2. + zc_T2**2.)*area_T2   
        MOI = MOI + I_T1 + I_T2 #need to multiply by rho*a or rhos*a to get mass MOI
        TotalArea = TotalArea + area_T1 + area_T2
    cx = cxSum/TotalArea/6.
    cz = czSum/TotalArea/6.
    wettedx = max(x) - min(x)
    xtil = (max(x) + min(x))/2.
    return TotalArea,cx,cz,MOI,wettedx,xtil

def overturning(w,t,p):
    """
    This function countains the model for block overturning. The function allows for 
    partially submerged and completely submerged blocks.
    """    
    from numpy import sqrt, pi, sin, cos, min
    th,dthdt = w
    a,b,c,alpha,rough ,uu, rhos, cd, cl, mu, eta, CmA, thetacr = p
    rho = 1032.0
    gra = 9.81
    lo = sqrt(0.25*(c**2.+b**2.))
    lW = lo*sin(thetacr-th) 
    Wt = rhos*gra*a*b*c
    IA = rhos*a*b*c*(c**2.+b**2.)/3.
    if eta < c: 
        lD = 0.5*eta

        if th<=thetacr + +pi/100.:
            xin = array([0.,-b,-b,0.,0.]) #vertex 0 included twice
            zin = array([0.,0.,c,c,0.]) #vertex 0 included twice
            xrot,zrot = rotatecoord(xin,zin,th)
            Asub,xc,zc,areaI,btil,xtil = polymoments(xrot,zrot,lo,eta)
        else:
            Asub = 0.
            xc = 0.
            zc = 0.
            areaI = 0.
            btil = 0.
            xtil = 0.

#         print('new:',th,Asub,xc,zc,areaI,btil,xtil)           

#cross-section option        Asub,xc,zc,areaI,btil,xtil = blockmoments(b,c,lo,eta,th,thetacr)
#polygon functions option         Asub,xc,zc,areaI,btil,xtil = shapemoments(b,c,lo,eta,th,thetacr)
        Vsub = a*Asub
        lB = -xc
        lL = -xtil
        IM = CmA*rho*a*areaI
    elif (lo*cos(thetacr-th) > 0.5*eta):
        lD = 0.5*eta

        if th<=thetacr + +pi/100.:
            xin = array([0.,-b,-b,0.,0.]) #vertex 0 included twice
            zin = array([0.,0.,c,c,0.]) #vertex 0 included twice
            xrot,zrot = rotatecoord(xin,zin,th)
            Asub,xc,zc,areaI,btil,xtil = polymoments(xrot,zrot,lo,eta)
        else:
            Asub = 0.
            xc = 0.
            zc = 0.
            areaI = 0.
            btil = 0.
            xtil = 0.

#         print('new:',th,Asub,xc,zc,areaI,btil,xtil)            

#cross-section option        Asub,xc,zc,areaI,btil,xtil = blockmoments(b,c,lo,eta,th,thetacr)
#polygon functionsoption         Asub,xc,zc,areaI,btil,xtil = shapemoments(b,c,lo,eta,th,thetacr)
        Vsub = a*Asub
        lB = -xc
        lL = -xtil
        IM = CmA*rho*a*areaI
    else: 
        lD = lo*cos(thetacr-th)
        lL = lW
        btil = b*cos(th) + c*sin(th)
        Vsub = a*b*c
        lB = lW
        IM = CmA*rho*a*b*c*(c**2.+b**2.)/3.
    FB = rho*gra*Vsub
    FL = 0.5*cl*rho*(uu**2.)*a*btil
    FD = 0.5*cd*rho*(uu**2.)*a*(2*lD)
    d2thdt2 = ( -Wt*lW + FB*lB + FL*lL + FD*lD ) / (IA + IM)
    sumM = -Wt*lW + FB*lB + FL*lL + FD*lD  #sum of moments	
    if sumM > 0.:
        f = [dthdt,d2thdt2]
    else:
        f = [0.,0.]
    return f
####------Robert------------------------------------------------


def rotate_boulder(x0,y0,angle_ll):
#    from numpy import append, array,zeros,cos,sin
#    import numpy as np
    # print(x0)
    # print(y0)
    # print(angle_ll)
    xP = zeros(len(x0))
    yP = zeros(len(x0))
    angle_l = angle_ll
    for i in range(len(x0)):
        xP[i] = x0[i]*cos(angle_l) - y0[i] * sin(angle_l)
        yP[i] = x0[i]*sin(angle_l) + y0[i] * cos(angle_l)  
    return xP, yP



def center_of_mass_xy(newbx_l,newby_l):
    """
    Calculating the (x,y) coordinates of the center of mass (COM) with equations
    proposed by Paul Bourke in 1988.
    """
#-------------------Area of polygon-----------------------------------
#     dummy = 0.0
#     for i in range(0,len(newbx_l)-1):
#         dummy = dummy + newbx_l[i] * newby_l[i+1] - newbx_l[i+1] * newby_l[i]
#     b_area_l = 0.5*dummy
    b_area_l = find_area(newbx_l,newby_l)
#------------------------------------------------------------------------------

#-------------------x coordinate of COM----------------------------------------
    dummy = 0.0
    for i in range(0,len(newbx_l)-1):
        dummy += (newbx_l[i] +
                newbx_l[i+1])*(newbx_l[i]*newby_l[i+1]-newbx_l[i+1]*newby_l[i])
    com_f_x_l = dummy/(6.0 * b_area_l)
#-------------------------------------------------------------------------------

#------------------y coordinate of COM-----------------------------------------
    dummy = 0.0
    for i in range(0,len(newbx_l)-1):
        dummy += (newby_l[i] +
                newby_l[i+1])*(newbx_l[i]*newby_l[i+1]-newbx_l[i+1]*newby_l[i])
    com_f_y_l = dummy/(6.0 * b_area_l)
#-------------------------------------------------------------------------------
    return com_f_x_l, com_f_y_l

def find_area(arr_x,arr_y):
    dummy = 0.0
    for i in range(0,len(arr_x)-1):
        dummy = dummy + arr_x[i] * arr_y[i+1] - arr_x[i+1] * arr_y[i]
    return 0.5*dummy


def find_new_cords(b_x_l,b_y_l,l_x_l,l_y_l):

##--------------Fully submerged-----------------------------------------------
    if l_y_l[0] > max(b_y_l):
        p_x = b_x_l
        p_y = b_y_l
    else:   
#---------Determine the faces for entry and exit-----------------------
        for i in range(0,len(b_x_l)-1):
##---------------Entry-------------------------------------------
            if l_y_l[0] >= b_y_l[i] and l_y_l[0] <= b_y_l[i+1]:
#                print  'Entry:',i, i+1
                entry_index = i
                break
        for i in range(0,len(b_x_l)):
##----------------Exit---------------------------------------------
            if l_y_l[0] <= b_y_l[i] and l_y_l[0] >= b_y_l[i+1]:
#                print  'Exit:', i, i+1
                exit_index = i
                break
#------------------------------------------------------------------

#--------------Determine coordinates for entry and exit points--------
##---------------------------Entry------------------------------------
        b_m_entry = (b_y_l[entry_index+1] - b_y_l[entry_index]) / (b_x_l[entry_index+1] - b_x_l[entry_index]) 
        b_b_entry = b_y_l[entry_index] - b_m_entry * b_x_l[entry_index]
        l_m_entry = (l_y_l[1] - l_y_l[0]) / (l_x_l[1] - l_x_l[0]) 
        l_b_entry = l_y_l[0] - l_m_entry * l_x_l[0]
        entry_intersect_x = (b_b_entry - l_b_entry) / (l_m_entry-b_m_entry)
        entry_intersect_y = b_m_entry * entry_intersect_x + b_b_entry
##---------------------------Exit--------------------------------------
        b_m_exit = (b_y_l[exit_index+1] - b_y_l[exit_index]) / (b_x_l[exit_index+1] - b_x_l[exit_index]) 
        b_b_exit = b_y_l[exit_index] - b_m_exit * b_x_l[exit_index]
        l_m_exit = (l_y_l[1] - l_y_l[0]) / (l_x_l[1] - l_x_l[0]) 
        l_b_exit = l_y_l[0] - l_m_exit * l_x_l[0]
        exit_intersect_x = (b_b_exit - l_b_exit) / (l_m_exit-b_m_exit)
        exit_intersect_y = b_m_exit * exit_intersect_x + b_b_exit
#------------------------------------------------------------------------

#-----------------------Make polygon-------------------------------------
        p_x = []
        p_y = []
        for i in range(0,entry_index+1):
            p_x.append(b_x_l[i])
            p_y.append(b_y_l[i])
        p_x.append(entry_intersect_x)
        p_y.append(entry_intersect_y)
        p_x.append(exit_intersect_x)
        p_y.append(exit_intersect_y)
        for i in range(exit_index+1,len(b_x_l)):
            p_x.append(b_x_l[i])
            p_y.append(b_y_l[i])
        p_x = array(p_x)
        p_y = array(p_y)
#----------------------------------------------------------------------------    
    return p_x,p_y



def run_depth_ode(a,b,c,froude,eta_v,rhos,cd,cl,mu,CmA,slope_l,delta_l,stoptime,nt):
    """ 
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers. To solve the equations, odespy is emplpoyed.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    import ode
    #rhos = float(rhos)
    roughness = 0.0 #JI set for Towers
#8    stoptime =10. #JI see Weiss & Diplas
 #   nt=1000
    t=linspace(0,stoptime,nt)
    dt = t[1]-t[0]
    alpha = arctan((0.5*b)/(0.5*c-delta_l))
    alpha_s = slope_l * pi /180.
    if delta_l != 0.5*c:
        ini_angle = 0.5* pi - arctan((0.5*b)/(0.5*c-delta_l))-alpha_s #JI 0.0 * pi/180.
    if delta_l == 0.5*c:
         ini_angle = 0.0 -alpha
    if delta_l> 0.5 *c:
        ini_angle = -(0.5*pi + arctan((0.5*b)/(0.5*c-delta_l)) + alpha_s)
    global thetacr
    thetacr = 0.5*pi - alpha_s#JI 0.5*pi - arctan(c/b) - alpha
    global p
    p = [a,b,c,alpha,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
#----------------Rotation---------------------------------------------------------
    w0_rot = [ini_angle,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    ti_rot, theta_rot = ode.euler(dis_2ndmoment_ode,w0_rot,[0,stoptime],dt)
# #---------------------------------------------------------------------------------
# 
# #----------------Sliding----------------------------------------------------------
    if delta_l== 0.0:
        w0_sli = [0.0,0.0]
        abserr = 1.0e-8
        relerr = 1.0e-6
        ti_sli, theta_sli = ode.euler(sliding_ode,w0_sli,[0,stoptime],dt)
        theta_rot = array(theta_rot)
        theta_sli = array(theta_sli)
        #print(shape(theta))
        return theta_rot[:,0],theta_sli[:,0]
    else:
        theta_rot = array(theta_rot)
        return theta_rot[:,0]

def run_depth_ode_1(a,b,c,uu_v,eta_v,rhos,cd,cl,mu,CmA,slope_l,delta_l,stoptime,nt):
    """ 
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers. To solve the equations, odespy is emplpoyed.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    
    import ode
    bx = [0.0,b,b,0.0,0.0]
    by = [-delta_l,-delta_l,c-delta_l,c-delta_l,-delta_l]
    t=linspace(0,stoptime,nt)
    dt = t[1]-t[0]
    sl_x= [-3.0,0.0,0.0,b+1]
    sl_y =[0.0,0.0,-delta_l,-delta_l]
    global p
    alpha_s = -slope_l * pi /180.
    bx,by=rotate_boulder(bx,by,alpha_s)
    rsl_x,rsl_y=rotate_boulder(sl_x,sl_y,alpha_s)
    com_obx,com_oby = center_of_mass_xy(bx,by)
    if com_oby >= 0.0:
        alpha = arctan((com_obx)/(com_oby))
        ini_angle = arctan((com_oby)/(com_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,bx,by,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
#     if com_oby == 0.0:
#         alpha = 0.0
#         ini_angle =  1.5* pi - arctan((com_oby)/(com_obx))
#         p = [alpha_s,bx,by,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
    if com_oby<0.0:
        abserr = 1.0e-8
        relerr = 1.0e-6
        nbx,nby = rotate_boulder(bx,by,-alpha_s)
        scom_obx,scom_oby = center_of_mass_xy(bx,by)
        w0_sli = [scom_oby,0.0]
        sini_angle = arctan((scom_oby)/(scom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,sini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
        ti_vsli, theta_vsli = ode.euler(vert_sliding_ode,w0_sli,[0,2*stoptime],dt)
        theta_vsli = asarray(theta_vsli)
        sl_index = 0
#        print(shape(theta_vsli))
        for i_sl in range(len(theta_vsli[0,:])):
#            print(i_sl,theta_vsli[0,i_sl])
            if theta_vsli[0,i_sl] > 0.0:
                sl_index = i_sl
                break
#        print('\t',sl_index,max(theta_vsli[0,:]),mu)
        for i in range(len(nbx)):
            nby[i] = nby[i]+ abs(theta_vsli[0,0]-theta_vsli[0,sl_index])
#        print(nby)
#        print(by)
        bx,by = rotate_boulder(nbx,nby,alpha_s)
        rcom_obx,rcom_oby = center_of_mass_xy(bx,by)
      #  print(scom_obx,scom_oby,rcom_obx,rcom_oby)
        com_obx = rcom_obx
        com_oby = rcom_oby
        ini_angle = arctan((rcom_oby)/(rcom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]

#----------------Rotation---------------------------------------------------------
  #  print(ini_angle*180./3.14)
    w0_rot = [ini_angle,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    global sumF
    sumF = [] 
    ti_rot, theta_rot = ode.euler(dis_2ndmoment_ode_dbk,w0_rot,[0,stoptime],dt)
    theta_rot = array(theta_rot)
###----sliding--------------------
    if delta_l== 0.0:
        scom_obx,scom_oby = center_of_mass_xy(bx,by)
        w0_sli = [scom_obx,0.0]
        abserr = 1.0e-8
        relerr = 1.0e-6
        ti_sli, theta_sli = ode.euler(sliding_ode,w0_sli,[0,stoptime],dt)
        theta_rot = array(theta_rot)
        theta_sli = array(theta_sli)
        #print(shape(theta))
        return theta_rot[0,:],theta_sli[0,:]
    else:
        theta_rot = array(theta_rot)
        dummy_local = theta_rot[0,:]
        return dummy_local

def run_depth_ode_1i(a,b,c,uu_v,eta_v,rhos,cd,cl,mu,CmA,slope_l,delta_l,stoptime,nt):
    """ 
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers. To solve the equations, odespy is emplpoyed.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    
    import ode
    bx = [0.0,b,b,0.0,0.0]
    by = [-delta_l,-delta_l,c-delta_l,c-delta_l,-delta_l]
    t=linspace(0,stoptime,nt)
    dt = t[1]-t[0]
    sl_x= [-3.0,0.0,0.0,b+1]
    sl_y =[0.0,0.0,-delta_l,-delta_l]
    global p
    alpha_s = -slope_l * pi /180.
    bx,by=rotate_boulder(bx,by,-alpha_s)
    rsl_x,rsl_y=rotate_boulder(sl_x,sl_y,-alpha_s)
    com_obx,com_oby = center_of_mass_xy(bx,by)
    if com_oby > 0.0:
        alpha = arctan((com_obx)/(com_oby))
        ini_angle = 0.5* pi - arctan((com_oby)/(com_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,bx,by,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
    if com_oby == 0.0:
        alpha = 0.0
        ini_angle = 0.0 +alpha_s
        p = [alpha_s,bx,by,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
    if com_oby<0.0:
        abserr = 1.0e-8
        relerr = 1.0e-6
        nbx,nby = rotate_boulder(bx,by,alpha_s)
        scom_obx,scom_oby = center_of_mass_xy(bx,by)
        w0_sli = [scom_oby,0.0]
        sini_angle = 0.5* pi - arctan((scom_oby)/(scom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,sini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]
        ti_vsli, theta_vsli = ode.euler(vert_sliding_ode,w0_sli,[0,2*stoptime],dt)
        theta_vsli = asarray(theta_vsli)
        sl_index = 0
#        print(shape(theta_vsli))
        for i_sl in range(len(theta_vsli[0,:])):
#            print(i_sl,theta_vsli[0,i_sl])
            if theta_vsli[0,i_sl] > 0.0:
                sl_index = i_sl
                break
#        print('\t',sl_index,max(theta_vsli[0,:]),mu)
        for i in range(len(nbx)):
            nby[i] = nby[i]+ abs(theta_vsli[0,0]-theta_vsli[0,sl_index])
#        print(nby)
#        print(by)
        bx,by = rotate_boulder(nbx,nby,-alpha_s)
        rcom_obx,rcom_oby = center_of_mass_xy(bx,by)
      #  print(scom_obx,scom_oby,rcom_obx,rcom_oby)
        com_obx = rcom_obx
        com_oby = rcom_oby
        ini_angle = arctan((rcom_oby)/(rcom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,ini_angle,uu_v,rhos,cd,cl,mu,eta_v,CmA]

#----------------Rotation---------------------------------------------------------
  #  print(ini_angle*180./3.14)
    w0_rot = [ini_angle,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    global sumF
    sumF = [] 
    ti_rot, theta_rot = ode.euler(dis_2ndmoment_ode_dbk,w0_rot,[0,stoptime],dt)
    theta_rot = array(theta_rot)
###----sliding--------------------
    if delta_l== 0.0:
        scom_obx,scom_oby = center_of_mass_xy(bx,by)
        w0_sli = [scom_obx,0.0]
        abserr = 1.0e-8
        relerr = 1.0e-6
        ti_sli, theta_sli = ode.euler(sliding_ode,w0_sli,[0,stoptime],dt)
        theta_rot = array(theta_rot)
        theta_sli = array(theta_sli)
        #print(shape(theta))
        return theta_rot[0,:],theta_sli[0,:]
    else:
        theta_rot = array(theta_rot)
        dummy_local = theta_rot[0,:]
        return dummy_local




def run_depth_ode_dbk(a,b,c,froude,eta_v,rhos,cd,cl,mu,CmA,slope_l,delta_l,stoptime,nt):
    """ 
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers. To solve the equations, odespy is emplpoyed.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    
    import ode
    fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(15,15)) 
    #rhos = float(rhos)
    # roughness = 0.0 #JI set for Towers
#8    stoptime =10. #JI see Weiss & Diplas
 #   nt=1000
    bx = [0.0,b,b,0.0,0.0]
    by = [-delta_l,-delta_l,c-delta_l,c-delta_l,-delta_l]
    t=linspace(0,stoptime,nt)
    dt = t[1]-t[0]
    sl_x= [-3.0,0.0,0.0,b+1]
    sl_y =[0.0,0.0,-delta_l,-delta_l]
    global p
#    plt.plot(bx,by)
    print(delta_l)
    alpha_s = slope_l * pi /180.
    bx,by=rotate_boulder(bx,by,-alpha_s)
    rsl_x,rsl_y=rotate_boulder(sl_x,sl_y,-alpha_s)
    com_obx,com_oby = center_of_mass_xy(bx,by)
    ax1.plot(bx,by)
    ax1.plot(com_obx,com_oby,'ko')
    ax1.plot(rsl_x,rsl_y,'k')
    ax1.axhline(0.0)
    ax1.set_ylim(-10,10)
    ax1.set_xlim(-10,10)
    if com_oby > 0.0:
        alpha = arctan((com_obx)/(com_oby))
        ini_angle = 0.5* pi - arctan((com_oby)/(com_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,bx,by,a,b,c,ini_angle,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
    if com_oby == 0.0:
        alpha = 0.0
        ini_angle = 0.0 +alpha_s
        p = [alpha_s,bx,by,a,b,c,ini_angle,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
    if com_oby<0.0:

        abserr = 1.0e-8
        relerr = 1.0e-6
        nbx,nby = rotate_boulder(bx,by,alpha_s)
        scom_obx,scom_oby = center_of_mass_xy(bx,by)
        w0_sli = [scom_oby,0.0]
        sini_angle = 0.5* pi - arctan((scom_oby)/(scom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,sini_angle,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
        ti_vsli, theta_vsli = ode.euler(vert_sliding_ode,w0_sli,[0,stoptime],dt)
        theta_vsli = asarray(theta_vsli)
        sl_index = 0
        for i_sl in range(len(theta_vsli[:,0])):
     #       print(i_sl,theta_vsli[i_sl,0])
            if theta_vsli[i_sl,0] > 0.0:
                sl_index = i_sl
                break
    #    print(scom_oby,sl_index,theta_vsli[0,0])        
        for i in range(len(nbx)):
            nby[i] = nby[i]+ abs(theta_vsli[0,0]-theta_vsli[sl_index,0])
        bx,by = rotate_boulder(nbx,nby,-alpha_s)
        rcom_obx,rcom_oby = center_of_mass_xy(bx,by)
        ax1.plot(rcom_obx,rcom_oby,'ro')
        com_obx = rcom_obx
        com_oby = rcom_oby
     #   print('shape results',shape(theta_vsli))
        ini_angle = arctan((rcom_oby)/(rcom_obx)) #JI 0.0 * pi/180.
        p = [alpha_s,nbx,nby,a,b,c,ini_angle,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
    ax1.plot(bx,by,'k:')
#         theta_vsli = asarray(theta_vsli)
#        print('test:',theta_vsli[-1,0])
#         if theta_vsli[-1,0] >= 0.0:
#                 for i in range(len(bx)):
#                     bx[i] = bx[i]
#                     by[i] = by[i]+ 0.5*c
#                bx = [0.0,b,b,0.0,0.0]
#                by = [0.5*c,0.5*c,c+0.5*c,c+0.5*c,0.5*c]
#                 p = [bx,by,a,b,c,alpha,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA]
#                 ini_angle = 0.0 -alpha_s
#    return ti_vsli,theta_vsli

#    global thetacr
#    thetacr = 0.5*pi - ini_angle#JI 0.5*pi - arctan(c/b) - alpha
    # 
    # print('initial angle:',ini_angle*180./pi)
    # print('Critical Angle:',thetacr*180./pi)

#----------------Rotation---------------------------------------------------------
    w0_rot = [ini_angle,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    global sumF
    sumF = [] 
    ti_rot, theta_rot = ode.euler(dis_2ndmoment_ode_dbk,w0_rot,[0,stoptime],dt)
    theta_rot = array(theta_rot)
    print(theta_rot[0,0],theta_rot[-1,0])
    index_dis = 0
    for i in range(len(theta_rot[:,0])):
        if theta_rot[i,0]> 0.5 * pi  -  arctan(com_oby/com_obx):
            index_dis=i
            break
    print(index_dis)
    print(index_dis,theta_rot[0,0],arctan(com_oby/com_obx))
    rbx,rby=rotate_boulder(bx,by,theta_rot[index_dis,0])
    print('min,max:',min(theta_rot[:,0]),max(theta_rot[:,0]))
    print('iniangle',0.5* pi - arctan((rcom_oby)/(rcom_obx)))
    com_dx,com_dy = center_of_mass_xy(rbx,rby)
    ax1.plot([0.0,com_dx],[0.0,com_dy],'ro')
    ax1.plot(rbx,rby,'r:')
    ax1.axhline(0.0)
    ax1.axvline(0.0)
#    print(shape(theta_rot))
#    print(theta_rot[-1,2],theta_rot[-1,0],theta_rot[-1,1])
#---------------------------------------------------------------------------------
# 
# #----------------Sliding----------------------------------------------------------
#     w0_sli = [0.0,0.0]
#     abserr = 1.0e-8
#     relerr = 1.0e-6
#     ti_sli, theta_sli = ode.euler(sliding_ode,w0_sli,[0,stoptime],dt)
# #----------------------------------------------------------------------------------
#     theta_rot = array(theta_rot)
#     theta_sli = array(theta_sli)
    # nbx,nby=rotate_boulder(bx,by,theta_rot[10,0]-ini_angle)
    # com_rox,com_roy = center_of_mass_xy(nbx,nby)
    # ax1.plot(nbx,nby)
    # ax1.plot(nbx[2],nby[2],'go')
    # ax1.plot(nbx[1],nby[1],'bo')
    # ax1.plot(com_rox,com_roy,'ro')
    # ax1.set_aspect('equal')
    # fig, ax2 = plt.subplots(nrows=1,ncols=1,figsize=(15,15)) 
    # ax2.plot(ti_rot,theta_rot[:,0],'k-')
    # print('y center of boulder:',com_oby)
    plt.show()
#    return ti_rot,theta_rot
    return 1,1

#


def vert_sliding_ode(t,w):
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1 = w
    slope_l,bx,by,aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA = p
    rho = 1032.0
    gra = 9.81
#    bx = [0.0,bb,bb,0.0,0.0]
#    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
# #----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(bb/cc)
    angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
# #------------------------------------------------------------
# 
#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)


    ld_v1 = com_a_y
    ll_h1 = com_a_x
    lb_h1 = com_a_x
    lw_h1 = com_f_x
    drag_c1 = float(abs(ad_y_min-ad_y_max))
    lift_b1 = float(abs(ad_x_min-ad_x_max))
    ld_v = ld_v1
    ll_h=ll_h1
    lb_h=lb_h1
    lw_h =lw_h1
    drag_c = drag_c1
    lift_b=lift_b1
    ll_f = (bb**2.0 + cc**2.0)
    if eta < 0.6*ad_y_max:
        ll_adj = (eta**2.0 + (vol_adj/eta)**2.0)
    else:
        ll_adj = ll_f
    m = ((1./3.*area_f*aa *rhos + CmA*rho* vol_adj))**(-1.0)
    fw = rhos * gra * vol_full
    fb = rho * gra * vol_adj
    cl = 1.0
#    fd = 0.5 * aa * drag_c * cd * rho * (uu-y1)**2.0
#    print(cl)
    fl=0.5*lift_b*aa*0.1*rho*(uu-y1)**2.0 #JI, not currently using lift, TODO if lift added: confirm bc is characteristic area and not ab
#    print(fl)\
#    print(mu)
    fsum2=fl +   (fb - fw)
    #print(fsum2,fl)
    term1 = fsum2 * m
    t1=t2=0.0
    if fl + (fb -fw) >0.0:
        t1 = y1
        t2 = term1
    return t1,t2    

    
def sliding_ode(t,w):
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1 = w
    slope_l,bx,by,aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA = p
    rho = 1032.0
    gra = 9.81
    # bx = [0.0,bb,bb,0.0,0.0]
    # by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
# #----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(bb/cc)
    angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
# #------------------------------------------------------------
# 
#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)


    ld_v1 = com_a_y
    ll_h1 = com_a_x
    lb_h1 = com_a_x
    lw_h1 = com_f_x
    drag_c1 = float(abs(ad_y_min-ad_y_max))
    lift_b1 = float(abs(ad_x_min-ad_x_max))
    ld_v = ld_v1
    ll_h=ll_h1
    lb_h=lb_h1
    lw_h =lw_h1
    drag_c = drag_c1
    lift_b=lift_b1
    ll_f = (bb**2.0 + cc**2.0)
    if eta < 0.6*ad_y_max:
        ll_adj = (eta**2.0 + (vol_adj/eta)**2.0)
    else:
        ll_adj = ll_f
    m = ((1./3.*area_f*aa *rhos + CmA*rho* vol_adj))**(-1.0)
    fw = rhos * gra * vol_full
    fb = rho * gra * vol_adj
    fd = 0.5 * aa * drag_c * cd * rho * (uu-y1)**2.0
#    fl=0.5*lift_b*aa*cl*rho*uu**2.0 #JI, not currently using lift, TODO if lift added: confirm bc is characteristic area and not ab
    fsum2=fd*ld_v + mu * (fb*lb_h - fw*lw_h)*cos(slope_l)
    term1 = fsum2 * m
    t1=t2=0.0
    if fd*ld_v + mu *(fb*lb_h -fw*lb_h)*cos(slope_l) >0.0:
        t1 = y1
        t2 = term1
    return t1,t2    
    
def dis_2ndmoment_ode(t,w):
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1 = w
    slope_l,bx,by,aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA = p
    #print(rhos,eta,CmA)
    rho = 1032.0
    gra = 9.81
#    bx = [0.0,bb,bb,0.0,0.0]
#    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
#----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(cc/bb)
    if angle == 0.0: angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
#------------------------------------------------------------

#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
    print(com_f_x,com_f_y)
    mi_f = mi_polygon(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
    mi_a = mi_polygon(adjbx,adjby)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)

    dummy_x, dummy_y =rotate_boulder(adjbx,adjby,slope_l)
    dc_x,dc_y = center_of_mass_xy(dummy_x,dummy_y)
    height = 0.5*(max(dummy_y))
    length = max(dummy_x)-0.5*(max(dummy_x)-min(dummy_x))
    ld_v = height
    ll_h = length
    lb_h = com_a_x
    lw_h = com_f_x
    drag_c = float(abs(ad_y_min-ad_y_max))
    lift_b = float(abs(ad_x_min-ad_x_max))
#     ll_f = (bb**2.0 + cc**2.0)
#     if eta < 0.6*ad_y_max:
#         ll_adj = (eta**2.0 + (vol_adj/eta)**2.0)
#     else:
#         ll_adj = ll_f
#     m= ((1./3. * ll_f*area_f*aa *rhos + ll_adj *CmA*rho* vol_adj*aa))**(-1.0)
    m = (mi_f * rhos + CmA * mi_a * rho)**(-1.0)
    fw=rhos*gra*vol_full
    fb=rho*gra*vol_adj
    fd=0.5*aa*drag_c*cd*rho*(uu-cos(y1))**2.0
    fl=0.5*lift_b*aa*cl*rho*(uu-cos(y1))**2.0 #JI, not currently using lift, TODO if lift added: confirm bc is characteristic area and not ab
    fsum2=fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h
    term1 = fsum2 * m
    t1=t2=0.0
    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lb_h >0.0:
        t1 = y1
        t2 =term1
    else:
        t1 = 0.0
        t2 = 0.0
#    if x1 < 0.5 * pi:
    return t1,t2

def dis_2ndmoment_ode_dbk(t,w):
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1= w
    slope_l,bx,by,aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA = p
    #print(rhos,eta,CmA)
    rho = 1032.0
    gra = 9.81
#    bx = [0.0,bb,bb,0.0,0.0]
#    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
#----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - alpha
    if angle == 0.0: angle=-0.00000000000001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
#------------------------------------------------------------

#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
#    print(com_f_x,com_f_y)
    mi_f = mi_polygon(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
    mi_a = mi_polygon(adjbx,adjby)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)
    
    dummy_x, dummy_y =rotate_boulder(adjbx,adjby,slope_l)
    dc_x,dc_y = center_of_mass_xy(dummy_x,dummy_y)
    height = 0.5*(max(dummy_y))
    length = max(dummy_x)-0.5*(max(dummy_x)-min(dummy_x))
    ld_v = height
    ll_h = length
    lb_h = com_a_x
    lw_h = com_f_x
    drag_c = float(abs(ad_y_min-ad_y_max))
    lift_b = float(abs(ad_x_min-ad_x_max))
    m = (mi_f * rhos + CmA * mi_a * rho)**(-1.0)
    fw=rhos*gra*vol_full
    fb=rho*gra*vol_adj
    fd=0.5*aa*drag_c*cd*rho*(uu-cos(y1))**2.0
    fl=0.5*lift_b*aa*cl*rho*(uu-cos(y1))**2.0 
    term3 = lb_h
    fsum2=fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h
    term1 = fsum2 * m
    t1=t2=0.0
    sumF.append(fb*lb_h - fw*lw_h)
    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h >0.0:
        t1 = y1
        t2 =term1
    else:
        t1 = 0.0
        t2 = 0.0
#    if x1 < 0.5 * pi:
    return t1,t2

def mi_polygon(arr_x,arr_y):
    """
    Function to calculate the moment of intertial based on the polygon and 
    perpendicular axis theorem: J_z = I_x + I_y
    """
    nx = len(arr_x)
    
    dummy = 0
    for i in range(nx-1):
        term1 = arr_x[i]*arr_y[i+1] - arr_x[i+1]*arr_y[i]
        term2 = arr_y[i]**2. + arr_y[i]*arr_y[i+1] + arr_y[i+1]**2.
        dummy += term2 * term1
    ix = 1./12. * dummy
    dummy = 0
    for i in range(nx-1):
        term1 = arr_x[i]*arr_y[i+1] - arr_x[i+1]*arr_y[i]
        term2 = arr_x[i]**2. + arr_x[i]*arr_x[i+1] + arr_x[i+1]**2.
        dummy += term2 * term1
    iy = 1./12. * dummy
    com_l_x,com_l_y = center_of_mass_xy(arr_x,arr_y)
    ll = sqrt(com_l_x**2.0 + com_l_y**2.0)
    area = find_area(arr_x,arr_y)
    return (ix+iy) + area*ll**2.0
    
#    dummy = 0.0
#    for i in range(0,nx-1):
#        term1 = arr_x[i]**2.0 + arr_x[i] * arr_x[i+1] + arr_x[i+1]**2.0
#        term2 = arr_x[i] * arr_y[i+1] - arr_x[i+1] * arr_y[i]
#        dummy += term1 * term2
#    return 1./12. * dummy




def random_floats(low, high, size):
#    from random import uniform
#    from numpy import zeros,array
    # returns an array of random numbers of size 'size' between 'low' and 'hihg'.
    output = zeros(size) 
    for i in range(size):
        output[i] = uniform(low,high)
    output=array(output)
    return output

#     # returns an array of random numbers of size 'size' between 'low' and 'hihg'.
#     return [random.uniform(low, high) for _ in xrange(size)]


def sliding(w,t,p):
    """
    This function countains the model for block sliding. The function allows for 
    partially submerged and completely submerged blocks.
    """
    from numpy import min
    rho = 1032.0
    gra = 9.81
    x,dxdt = w
    a,b,c,alpha,rough ,uu, rhos, cd, cl, mu, eta, CmA, thetacr = p
    Cmx = 1.0
    mass = rhos*a*b*c
    addedmass = Cmx*rho*a*b*min([c,eta])
    W = gra*mass
    FB = rho*gra*a*b*min([c,eta])
    FL = 0.5*cl*rho*((uu - dxdt)**2.)*a*b
    FD = 0.5*cd*rho*((uu - dxdt)**2.)*a*min([c,eta])   
    d2xdt2 = (FD - mu*(W-FB-FL))/(mass+addedmass)
    FL0 = 0.5*cl*rho*(uu**2.)*a*b
    FD0 = 0.5*cd*rho*(uu**2.)*a*min([c,eta])
    sumFx = (FD0 - mu*(W-FB-FL0))
    if sumFx > 0.:
        f = [dxdt,d2xdt2]
    else:
        f = [0.,0.]
    return f



##JI modified version for Caesarea tower
def run_depth(froude,eta_v,a,b,c,rhos,cd,cl,mu,CmA):
    """
    Runnning the Monte-Carlo type simulations of block overturning and sliding for 
    different flow depths and Froude numbers.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    from scipy.integrate import odeint
    x1 = 0.5
    y1 = 0.0
    roughness = 0.0
    stoptime =5.5
    nt=1000
    t=linspace(0,stoptime,nt)
    alpha = 0.
    thetacr = arctan(b/c)
    p = [a,b,c,alpha,roughness,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA,thetacr]    
    w0 = [0.0,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    #wsol = odeint(overturning, w0,t, args=(p,), atol = abserr, rtol = relerr)
    wsol = odeint(overturning, w0,t, args=(p,))
    #wsol = odeint(dislodgement, w0,t, args=(p,), atol = abserr, rtol = relerr)
#    wsol = odeint(dislodgement_rw, w0,t, args=(p,))
    theta = zeros(len(wsol[:,0]))
    theta[:] = wsol[:,0]
    wsolslide = odeint(sliding,w0,t,args=(p,))    
    xslide = zeros(len(wsolslide[:,0]))
    xslide[:] = wsolslide[:,0]
    if max(theta)> thetacr:
        tover = min(t[theta>thetacr])
#         return 1, xslide[t==tover], tover
        return 1, max(xslide), tover
    else:
        return 0, max(xslide), 0.

def create_data_randomblocks(fname,aaa,bbb_min,bbb_max,ccc_min,ccc_max,size_MC,n_pr,rhos_min,rhos_max,fd_min,fd_max,fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,cma_max,num_time,e_time):
    from numpy import linspace,sqrt,zeros
    from tqdm import tqdm_notebook
    from random import  uniform
    from multiprocessing import Pool
    n_coeff = 10
    total_size_MC = int(size_MC)
    Fr = random_floats(fr_min, fr_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
    bbb = random_floats(bbb_min,bbb_max, total_size_MC)
    ccc = random_floats(ccc_min,ccc_max, total_size_MC)
    print("Running  code...")
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global froude,flow_depth,cls,cds
    aaa_g = aaa
    print(len(rhos),rhos_min)
    bbb_g = bbb
    ccc_g = ccc
#     print(aaa_g,bbb_g,ccc_g)
    mu_g = mu
    CmA_g = CmA
    froude=[]
    froude = Fr
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    rhos_g = rhos
    cls = cl_arr
    cds = cd_arr
    global dummy
    dummy = []
    dummy2 = []
    dummy3 = []
    #dummy = zeros(total_size_MC)
    pool = Pool(processes = int(n_pr))
    inputs = range(total_size_MC)
    for res in tqdm_notebook(pool.imap(func_randomblocks, inputs),total=total_size_MC):
     #   m+=1
        #dummy[m]=res
        dummy.append(res[0])
        dummy2.append(res[1])
        dummy3.append(res[2])
#JI         dummy.append(res) #this was a merge error, see 3 lines up.
    pool.close()
    pool.terminate()
    fr_values1 = []
    fr_values2 = []
    eta_values1 = []
    eta_values2 = []
    rhos_values1 = []
    rhos_values2 = []
    CD_values1 = []
    CD_values2 = []
    CL_values1 = []
    CL_values2 = []
    CmA_values1 = []
    CmA_values2 = []
    mu_values1 = []
    mu_values2 = []
    xslide_values1 = []
    xslide_values2 = []
    tover_values1 = []
    tover_values2 = []
    b_values1 = []
    b_values2 = []
    c_values1 = []
    c_values2 = []
    for jj in xrange(total_size_MC):
        if dummy[jj] == 1:
            fr_values1.append(froude[jj])
            eta_values1.append(flow_depth[jj])
            rhos_values1.append(rhos[jj])
            CD_values1.append(cd_arr[jj])
            CL_values1.append(cl_arr[jj])
            CmA_values1.append(CmA[jj])
            mu_values1.append(mu[jj])
            xslide_values1.append(dummy2[jj])
            tover_values1.append(dummy3[jj])
            b_values1.append(bbb[jj])
            c_values1.append(ccc[jj])
        if dummy[jj] == 0:
            fr_values2.append(froude[jj])
            eta_values2.append(flow_depth[jj])
            rhos_values2.append(rhos[jj])
            CD_values2.append(cd_arr[jj])
            CL_values2.append(cl_arr[jj])
            CmA_values2.append(CmA[jj])
            mu_values2.append(mu[jj])
            xslide_values2.append(dummy2[jj])
            tover_values2.append(dummy3[jj])
            b_values2.append(bbb[jj])
            c_values2.append(ccc[jj])
    fname1 = "{t1}_{t2}.dat".format(t1=fname,t2="1")
    fname2 = "{t1}_{t2}.dat".format(t1=fname,t2="2")
    save_data(zip(fr_values1,eta_values1,rhos_values1,CD_values1,CL_values1,b_values1,c_values1,CmA_values1,mu_values1,xslide_values1,tover_values1),fname1)
    save_data(zip(fr_values2,eta_values2,rhos_values2,CD_values2,CL_values2,b_values2,c_values2,CmA_values2,mu_values2,xslide_values2,tover_values2),fname2)


def create_data(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,fd_max,fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,cma_max,num_time,e_time):
#OLD def create_data(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,fd_max,fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu,cma_min,cma_max):
    from numpy import linspace,sqrt,zeros
    from tqdm import tqdm_notebook
    from random import  uniform
    from multiprocessing import Pool
    n_coeff = 10
    total_size_MC = int(size_MC)
    Fr = random_floats(fr_min, fr_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
    print("Running  code...")
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global froude,flow_depth,cls,cds
    aaa_g = aaa
    print(len(rhos),rhos_min)
    bbb_g = bbb
    ccc_g = ccc
    print(aaa_g,bbb_g,ccc_g)
    mu_g = mu
    CmA_g = CmA
    froude=[]
    froude = Fr
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    rhos_g = rhos
    cls = cl_arr
    cds = cd_arr
    global dummy
    dummy = []
    dummy2 = []
    dummy3 = []
    #dummy = zeros(total_size_MC)
    pool = Pool(processes = int(n_pr))
    inputs = range(total_size_MC)
    for res in tqdm_notebook(pool.imap(func, inputs),total=total_size_MC):
     #   m+=1
        #dummy[m]=res
        dummy.append(res[0])
        dummy2.append(res[1])
        dummy3.append(res[2])
#JI         dummy.append(res) #this was a merge error, see 3 lines up.
    pool.close()
    pool.terminate()
    fr_values1 = []
    fr_values2 = []
    eta_values1 = []
    eta_values2 = []
    rhos_values1 = []
    rhos_values2 = []
    CD_values1 = []
    CD_values2 = []
    CL_values1 = []
    CL_values2 = []
    CmA_values1 = []
    CmA_values2 = []
    mu_values1 = []
    mu_values2 = []
    xslide_values1 = []
    xslide_values2 = []
    tover_values1 = []
    tover_values2 = []
    for jj in xrange(total_size_MC):
        if dummy[jj] == 1:
            fr_values1.append(froude[jj])
            eta_values1.append(flow_depth[jj])
            rhos_values1.append(rhos[jj])
            CD_values1.append(cd_arr[jj])
            CL_values1.append(cl_arr[jj])
            CmA_values1.append(CmA[jj])
            mu_values1.append(mu[jj])
            xslide_values1.append(dummy2[jj])
            tover_values1.append(dummy3[jj])
        if dummy[jj] == 0:
            fr_values2.append(froude[jj])
            eta_values2.append(flow_depth[jj])
            rhos_values2.append(rhos[jj])
            CD_values2.append(cd_arr[jj])
            CL_values2.append(cl_arr[jj])
            CmA_values2.append(CmA[jj])
            mu_values2.append(mu[jj])
            xslide_values2.append(dummy2[jj])
            tover_values2.append(dummy3[jj])
    fname1 = "{t1}_{t2}.dat".format(t1=fname,t2="1")
    fname2 = "{t1}_{t2}.dat".format(t1=fname,t2="2")
    save_data(zip(fr_values1,eta_values1,rhos_values1,CD_values1,CL_values1,CmA_values1,mu_values1,xslide_values1,tover_values1),fname1)
    save_data(zip(fr_values2,eta_values2,rhos_values2,CD_values2,CL_values2,CmA_values2,mu_values2,xslide_values2,tover_values2),fname2)


def func_randomblocks(imc):
    dummy1,xslide,tover=run_depth(froude[imc],flow_depth[imc],aaa_g,bbb_g[imc],ccc_g[imc],rhos_g[imc],cds[imc],cls[imc],mu_g[imc],CmA_g[imc])
    return dummy1,xslide,tover

def func(imc):
#    dummy1,xslide,tover=run_depth(froude[imc],flow_depth[imc],aaa_g,bbb_g,ccc_g,rhos_g[imc],cds[imc],cls[imc],mu_g,CmA_g,num_time_g,e_time_g)
    dummy1 = run_depth_ode_1(aaa_g,bbb_g,ccc_g,uu_g[imc],\
            flow_depth[imc],rhos_g[imc],cds[imc],cls[imc],\
            mu_g[imc],CmA_g[imc],slope_g,delta_g,e_time_g,num_time_g)

    return dummy1#,xslide,tover


def create_data_pp(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,\
    fd_max,fu_min,fu_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,\
    cma_max,slope_l,delta_l,num_time,e_time):
    from numpy import linspace,sqrt,zeros
    from tqdm import tqdm_notebook
    from random import  uniform
    from multiprocessing import Pool
    n_coeff = 10
    total_size_MC = int(size_MC)
    uu = random_floats(fu_min, fu_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
    print("Running  code...")
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global uu_g,flow_depth,cls,cds,num_time_g,e_time_g
    global slope_g,delta_g
    aaa_g = aaa
    print(len(rhos),rhos_min)
    bbb_g = bbb
    ccc_g = ccc
    print(aaa_g,bbb_g,ccc_g)
    mu_g = mu
    CmA_g = CmA
    uu_g=[]
    uu_g = uu
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    num_time_g = int(num_time)
    e_time_g = e_time
    rhos_g = rhos
    delta_g = delta_l
    slope_g = slope_l
    cls = cl_arr
    cds = cd_arr
    data_dis = zeros([int(num_time),total_size_MC+1])
    data_sli = zeros([int(num_time),total_size_MC+1])
    global dummy
    data_dis[:,0] = linspace(0.0,e_time,num_time)
    data_sli[:,0] = linspace(0.0,e_time,num_time)
    dummy = []
    pool = Pool(processes = int(n_pr))
    inputs = range(total_size_MC)
    for res in tqdm.tqdm(pool.imap(func, inputs),total=total_size_MC):
        #dummy = concatenate((dummy,res),axis=0)
     #   print(res)
        res = array(res)
        dummy.append(res)
        
    pool.close()
    pool.terminate()
    dummy = asarray(dummy)
    print(shape(data_dis))
    print('shape of dummy:',shape(dummy),len(shape(dummy)))
    if len(shape(dummy)) == 2:
        for i in range(total_size_MC):
     #      print(shape(data_dis),shape(dummy))
           data_dis[:,i+1] = array(dummy[i,:])
    #    print()
        fname1 = "{t1}_rotation.dat".format(t1=fname)
        fname3 = "{t1}_conditions.dat".format(t1=fname)
        save_md_data(fname1,data_dis)
        save_md_data(fname3,transpose([uu[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname3))
    if len(shape(dummy)) == 3:
        print(shape(dummy))
        print(shape(data_dis))
        print(shape(data_sli))
        print(shape(dummy))
        for i in range(total_size_MC):
           data_dis[:,i+1] = array(dummy[i,0,:])
           data_sli[:,i+1] = array(dummy[i,1,:])
        print()
        fname1 = "{t1}_rotation.dat".format(t1=fname)
        fname2 = "{t1}_sliding.dat".format(t1=fname)
        fname3 = "{t1}_conditions.dat".format(t1=fname)
        save_md_data(fname1,data_dis)
        save_md_data(fname2,data_sli)
        save_md_data(fname3,transpose([uu[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname2))
        print('\t {t1}'.format(t1=fname3))


def create_data_cl(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,\
    fd_max,fu_min,fu_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,\
    cma_max,slope_l,delta_l,num_time,e_time):
    from numpy import linspace,sqrt,zeros
    from tqdm import tqdm_notebook
    from random import  uniform
    from multiprocessing import Pool
    n_coeff = 10
    total_size_MC = int(size_MC)
    uu = random_floats(fu_min, fu_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
    print("Running  code...")
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global uu_g,flow_depth,cls,cds,num_time_g,e_time_g
    global slope_g,delta_g
    aaa_g = aaa
    print(len(rhos),rhos_min)
    bbb_g = bbb
    ccc_g = ccc
    print(aaa_g,bbb_g,ccc_g)
    mu_g = mu
    CmA_g = CmA
    uu_g=[]
    uu_g = uu
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    num_time_g = int(num_time)
    e_time_g = e_time
    rhos_g = rhos
    delta_g = delta_l
    slope_g = slope_l
    cls = cl_arr
    cds = cd_arr
    data_dis = zeros([int(num_time),total_size_MC+1])
    data_sli = zeros([int(num_time),total_size_MC+1])
    global dummy
    data_dis[:,0] = linspace(0.0,e_time,num_time)
    data_sli[:,0] = linspace(0.0,e_time,num_time)
    dummy = []
    pool = Pool(processes = int(n_pr))
    inputs = range(total_size_MC)
#     for res in tqdm.tqdm(pool.imap(func, inputs),total=total_size_MC):
#         #dummy = concatenate((dummy,res),axis=0)
#      #   print(res)
#         res = array(res)
#         dummy.append(res)
#         
    dummy = pool.map(func, inputs)
    pool.close()
    pool.terminate()
    dummy = asarray(dummy)
    print(shape(data_dis))
    print('shape of dummy:',shape(dummy),len(shape(dummy)))
    if len(shape(dummy)) == 2:
        for i in range(total_size_MC):
     #      print(shape(data_dis),shape(dummy))
           data_dis[:,i+1] = array(dummy[i,:])
    #    print()
        fname1 = "{t1}_rotation.dat".format(t1=fname)
        fname3 = "{t1}_conditions.dat".format(t1=fname)
        save_md_data(fname1,data_dis)
        save_md_data(fname3,transpose([uu[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname3))
    if len(shape(dummy)) == 3:
        print(shape(dummy))
        print(shape(data_dis))
        print(shape(data_sli))
        print(shape(dummy))
        for i in range(total_size_MC):
           data_dis[:,i+1] = array(dummy[i,0,:])
           data_sli[:,i+1] = array(dummy[i,1,:])
        print()
        fname1 = "{t1}_rotation.dat".format(t1=fname)
        fname2 = "{t1}_sliding.dat".format(t1=fname)
        fname3 = "{t1}_conditions.dat".format(t1=fname)
        save_md_data(fname1,data_dis)
        save_md_data(fname2,data_sli)
        save_md_data(fname3,transpose([uu[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname2))
        print('\t {t1}'.format(t1=fname3))



def create_data_arc(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,\
    fd_max,fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,\
    cma_max,num_time,e_time):
    from numpy import linspace,sqrt,zeros
    from random import  uniform
    from multiprocessing import Pool
    n_coeff = 10
    total_size_MC = int(size_MC)
    Fr = random_floats(fr_min, fr_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global froude,flow_depth,cls,cds,num_time_g,e_time_g
    aaa_g = aaa
    print(len(rhos),rhos_min)
    bbb_g = bbb
    ccc_g = ccc
    print(aaa_g,bbb_g,ccc_g)
    mu_g = mu
    CmA_g = CmA
    froude=[]
    froude = Fr
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    num_time_g = int(num_time)
    e_time_g = e_time
    rhos_g = rhos
    cls = cl_arr
    cds = cd_arr
    data_dis = zeros([int(num_time),total_size_MC+1])
    data_sli = zeros([int(num_time),total_size_MC+1])
    global dummy
    data_dis[:,0] = linspace(0.0,e_time,num_time)
    data_sli[:,0] = linspace(0.0,e_time,num_time)
    dummy = []
    time1 = time.time()
    pool = Pool(processes = int(n_pr))
    inputs = range(total_size_MC)
    dummy = pool.map(func, inputs)
        # for res in tqdm.tqdm(pool.imap(func, inputs),total=total_size_MC):
        # dummy.append(res)
    pool.close()
    pool.terminate()
    print('Runtime1:',time.time()-time1)
    dummy = asarray(dummy)
    for i in range(total_size_MC):
       data_dis[:,i+1] = array(dummy[i,0,:])
       data_sli[:,i+1] = array(dummy[i,1,:])
#     fname1="{t1}_rotation.dat".format(t1=fname)
#     fname2="{t1}_sliding.dat".format(t1=fname)
#     save_data(zip(froude[:],eta[:],dummy[:,0],dummy[:,1]),fname1)
#     save_data(zip(froude[:],eta[:],dummy[:,2],dummy[:,3]),fname2)
    print('Runtime2:',time.time()-time1)
    print('Saving data')
    fname1 = "{t1}_rotation.dat".format(t1=fname)
    fname2 = "{t1}_sliding.dat".format(t1=fname)
    fname3 = "{t1}_conditions.dat".format(t1=fname)
    save_md_data(fname1,data_dis)
    save_md_data(fname2,data_sli)
    save_md_data(fname3,transpose([froude[:],eta[:]]))
    print('Runtime3:',time.time()-time1)

def create_data_ipad(fname,aaa,bbb,ccc,size_MC,n_pr,rhos_min,rhos_max,fd_min,\
    fd_max,fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,cma_min,\
    cma_max,slope_l,delta_l,num_time,e_time):
    #import console
    n_coeff = 10
    total_size_MC = int(size_MC)
    Fr = random_floats(fr_min, fr_max, total_size_MC)
    cl_arr = random_floats(cl_min,cl_max,total_size_MC) 
    cd_arr = random_floats(cd_min,cd_max,total_size_MC) 
    eta = random_floats(fd_min, fd_max,total_size_MC)
    rhos = random_floats(rhos_min, rhos_max, total_size_MC)
    CmA = random_floats(cma_min,cma_max, total_size_MC)
    mu = random_floats(mu_min,mu_max, total_size_MC)
#    console.set_color(0.0,0.5,1.0)
    print("Running  code...")
 #   console.set_color()
    global aaa_g,bbb_g,ccc_g,rhos_min_g,rhos_max_g,mu_g,CmA_g,rhos_g
    global froude,flow_depth,cls,cds,slope_g,delta_g
    aaa_g = aaa
    bbb_g = bbb
    ccc_g = ccc
    mu_g = mu
    CmA_g = CmA
    froude=[]
    froude = Fr
    flow_depth = eta
    rhos_min_g = rhos_min
    rhos_max_g = rhos_max
    rhos_g = rhos
    cls = cl_arr
    cds = cd_arr
    delta_g = delta_l
    slope_g = slope_l
#    dummy = zeros(num_time,2)
    print(num_time,total_size_MC)
    data_dis = zeros([num_time,total_size_MC+1])
    data_sli = zeros([num_time,total_size_MC+1])
    # time1 = time.time()
    data_dis[:,0] = linspace(0.0,e_time,num_time)
    data_sli[:,0] = linspace(0.0,e_time,num_time)
#    console.set_color(1.0,0.0,0.0)    
    print('\tProcessing:')
    for i in range(total_size_MC):
        print(f"\r\t\t\t{i+1}/{total_size_MC}",end="")
        dummy=run_depth_ode_1(aaa_g,bbb_g,ccc_g,froude[i],flow_depth[i],rhos_g[i],cds[i],cls[i],mu_g[i],CmA_g[i],slope_g,delta_g,e_time,num_time)
        dummy = asarray(dummy)
        if delta_l == 0.0:
            data_dis[:,i+1] = array(dummy[0,:])
            data_sli[:,i+1] = array(dummy[1,:])
        else:    
            data_dis[:,i+1] = array(dummy[:])
   #     sys.stdout.flush()
   # print()
#     console.set_color()
#   dummy = array(dummy)
#     console.set_color(1.0,1.0,0.0)
#     print("\nRuntime (n = {t1}): {t2:04.1f}s.".format(t1=total_size_MC,t2= time.time()-time1))
#     print("\t\t\t".format(t2= time.time()-time1))
#    console.set_color()
    print()
    fname1 = "{t1}_rotation.dat".format(t1=fname)
    fname2 = "{t1}_sliding.dat".format(t1=fname)
    fname3 = "{t1}_conditions.dat".format(t1=fname)
    if delta_l == 0.0:
        save_md_data(fname1,data_dis)
        save_md_data(fname2,data_sli)
        save_md_data(fname3,transpose([froude[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname2))
        print('\t {t1}'.format(t1=fname3))

    else:    
        save_md_data(fname1,data_dis)
        save_md_data(fname3,transpose([froude[:],eta[:]]))
        print('Results saved in:') 
        print('\t {t1}'.format(t1=fname1))
        print('\t {t1}'.format(t1=fname3))
