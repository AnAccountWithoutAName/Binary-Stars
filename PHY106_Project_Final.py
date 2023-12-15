#************************************************************* MODULES IMPORTED ******************************************************
from numpy import arange,real
from cmath import pi,sin,cos,acos,sqrt,atan
from math import radians
from datetime import *
from time import *
import matplotlib.pyplot as pl
import math as m
from SimulatorFinale import Simulate
#**************************************************************** ORBIT CALCULATOR ***************************************************

def StarType(L):
    lambda_peak=float(input("Enter the peak wavelength of the Star. "))
    b=2.8981e-3
    T=b/(lambda_peak*5778)      #Peak Surface Temperature in terms of Solar Temperature
    print('Peak Surface Temperature of the star is : ',T,' (in terms of Solar Temperature)')
    R=sqrt(L)* (1/(T**2))       #Radius of the Star in terms of Solar Radius
    print('Radius of the star in terms of Solar Radius is: ',R.real)
    
    if R.real<0.1 and L<1:
        print('The star is a white dwarf.')
    elif R.real>=0.1 and R.real<1 and L>=1 and L<10**2:
        print('The star is a main sequence star.')
    elif R.real>=1 and R.real<10:
        print('Star is of similar magnitude as that of Sun.')
    elif R.real>=10 and R.real<=100 and L>=10**2 and L<10**3:
        print('The star is a Red Giant.')
    elif R.real>=10 and R.real<=1000 or R.real>1000 and L>=10**4 and L<=10**6:
        print('The star is a supergiant.')
        
    if T>=1.7307 and T<=6.9228:
        print('Star is Hot.')
    elif T>=0.4326 and T<=0.8653:
        print('Star is Cool.')
    if L>100:
        print('Star is Bright.')
    elif L>=0.01 and L<=100:
        print('Star is neither too bright not too dim.')
    elif L<0.01: 
        print('Star is Dim.')
    
    n=int(input('Enter the max limit for distance away from the star.'))
    D=[]
    for i in range(1,n+1):
        D.append(i)
    sigma=float(input('Enter the value of Stefan-Boltzmann Constant for the star.'))
    Temp=[]
    for i in range(1,n+1):
       tem=(L/(16*pi*sigma*i**2))**(1/4) 
       Temp.append(tem)
    print('TEMPERATURE = ', Temp)
    
    return D, Temp
    
    
def ImgPlot(image_file):
    
    import numpy as np
    from PIL import Image
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from astropy.visualization import astropy_mpl_style
    
    default='sirius.jpg'
    if image_file==default:
         ##############################################################################
         # Use `astropy.io.fits.info()` to display the structure of the file:
             
             fits.info('red.fits')
             
    
        # Generally the image information is located in the Primary HDU, also known
        # as extension 0. Here, we use `astropy.io.fits.getdata()` to read the image
        # data from this first extension using the keyword argument ``ext=0``:

             image_data = fits.getdata('red.fits', ext=0)
            
            # The data is now stored as a 2D numpy array. Print the dimensions using the
            # shape attribute:

             print(image_data.shape)
    
        # Display the image data:
        
             plt.figure()
             plt.imshow(image_data, cmap='gnuplot')
             plt.colorbar()
    else: 
       
        plt.style.use(astropy_mpl_style)
    
        # Load and display the original 3-color jpeg image:
            
        image = Image.open(image_file)
        xsize, ysize = image.size
        print(f"Image size: {ysize} x {xsize}")
        print(f"Image bands: {image.getbands()}")
        ax = plt.imshow(image)
            
        # Split the three channels (RGB) and get the data as Numpy arrays. The arrays
        # are flattened, so they are 1-dimensional:
                
        r, g, b = image.split()
        r_data = np.array(r.getdata()) # data is now an array of length ysize*xsize
        g_data = np.array(g.getdata())
        b_data = np.array(b.getdata())
        print(r_data.shape)
                
        # Reshape the image arrays to be 2-dimensional:
        
        r_data = r_data.reshape(ysize, xsize) # data is now a matrix (ysize, xsize)
        g_data = g_data.reshape(ysize, xsize)
        b_data = b_data.reshape(ysize, xsize)
        print(r_data.shape)
                    
        # Write out the channels as separate FITS images.
        # Add and visualize header info
                        
        red = fits.PrimaryHDU(data=r_data)
        red.header['LATOBS'] = "32:11:56" # add spurious header info
        red.header['LONGOBS'] = "110:56"
        red.writeto('newred.fits')

        green = fits.PrimaryHDU(data=g_data)
        green.header['LATOBS'] = "32:11:56"
        green.header['LONGOBS'] = "110:56"
        green.writeto('newgreen.fits')
    
        blue = fits.PrimaryHDU(data=b_data)
        blue.header['LATOBS'] = "32:11:56"
        blue.header['LONGOBS'] = "110:56"
        blue.writeto('newblue.fits')

        from pprint import pprint
        pprint(red.header)
    

        ##############################################################################
        # Use `astropy.io.fits.info()` to display the structure of the file:

        fits.info('newred.fits')
            
            
        # Generally the image information is located in the Primary HDU, also known
        # as extension 0. Here, we use `astropy.io.fits.getdata()` to read the image
        # data from this first extension using the keyword argument ``ext=0``:

        image_data = fits.getdata('newred.fits', ext=0)
                
        # The data is now stored as a 2D numpy array. Print the dimensions using the
        # shape attribute:

        print(image_data.shape)
    
        # Display the image data:
        
        plt.figure()
        plt.imshow(image_data, cmap='gnuplot')
        plt.colorbar()

    

def radial_velocity(pds,argperi_deg,m1,m2): 
     k1=((((m2*sin(argperi_deg))**3)*2*pi*4.3009e-3)/(((m1+m2)**2)*pds))**(1/3)
     k2=((((m1*sin(argperi_deg)**3)*2*pi*4.3009e-3)/(((m1+m2)**2)*pds)))**(1/3)
     omega=2*pi/pds
     V1=[]
     V2=[]
     for j in range(0,int(pds)):
         v1=k1*m.sin(omega*j+argperi_deg)
         v2=k2*m.sin(omega*j+argperi_deg)
         V1.append(v1)
         V2.append(v2)
     
     return V1,V2
     

    
    
    
def orbit(t,sys_dys,ang_arcsec,inc_deg,longAN_deg,e,T,Pds,argperi_deg,m1,m2,alpha,alpha1):
    
    #MEAN MOTION
    n=2*pi/Pds
    
    M_deg=n*(Time_period_of_obs)
    #MEAN ANOMALY
    M_deg=radians(M_deg)

    #ECCENTRIC ANOMALY
    e_last=M_deg+ e*sin(M_deg)+ (e**2/(2*M_deg))*sin(2*M_deg)
    
    for i in range(4):
        M0_deg=e_last-e*sin(e_last)
        e_latest_last=e_last +((M_deg-M0_deg)/(1- e*cos(e_last)))
        e_last=e_latest_last
    
    #TRUE ANOMALY
    TA=acos((cos(e_latest_last)-e)/(1-e*cos(e_latest_last)))
    
    #Seperation Distance in Reference Plane from provided Seperation Distance and Seperation Angle AU
    sep_dis_rp=ang_arcsec*sys_dys
    
    #Seperation Distance in Orbital Plane
    sep_dis_op=sep_dis_rp/(sqrt(((sin(argperi_deg+TA)*((((cos(complex(0,1))).real))**2)+(cos((argperi_deg+TA))**2))))).real
    
    #Seperation Distance in Reference Plane with respect to the Line of Nodes
    sep_dis_rp_y=sep_dis_op*sin(argperi_deg+TA)*((cos(complex(0,1))).real)       #Y-component
    sep_dis_rp_x=sep_dis_op*cos(argperi_deg+TA)                                  #X-component
    
    #Angle between the Line of Nodes and Primary Inclination
    angln_pi_deg=atan(sep_dis_rp_y/sep_dis_rp_x)
    
    #Position Angle in Image 
    quad=int(input("Enter the quadrant in which primary star(m1) Lies => "))
    if quad==1 or quad==4:
        pos_angle=longAN_deg+angln_pi_deg +180
    else:
        pos_angle=longAN_deg+angln_pi_deg
        
    #Sum of SemiMajor axes of m1 and m2
    a_total=sep_dis_op*(1+e*cos(TA))/(1-e**2)
    
    #SemiMajor Axis of m1
    a1=a_total/(1+(m2/m1))
    
    #SemiMajor Axis of m2
    a2=(a1/m2)*m1
    
    #Radial Velocity
    v1r=((((m2*sin(inc_deg))**3)*2*pi*4.3009e-3)/(((m1+m2)**2)*Pds))**(1/3)
    v2r=((((m1*sin(inc_deg)**3)*2*pi*4.3009e-3)/(((m1+m2)**2)*Pds)))**(1/3)
    
    #Luminosity of m1 and m2
    L1=((m1)**alpha)
    L2=((m2)**alpha1)
    
    return a1,a2,v1r,v2r,L1,L2,pos_angle,TA,e_latest_last,M_deg,n

#**************************************************************** INPUTS *************************************************************
print('/n *********************** INPUTS *********************** /n')
print('/n *********************** INPUTS *********************** /n')
t=float(input("Enter epoch of observation (in epoch timestamp)= ")) 
sys_dys=float(input("Enter system distance from earth (in parsec)= ")) 
ang_arcsec=float(input("Enter seperation angle of stars (in degrees) = "))
inc_deg=float(input("Enter inclination angle (in degrees) = "))
longAN_deg=float(input("Enter longitude of ascending node (in degrees) = "))
e=float(input("Enter eccentricity of the orbit = "))
T=float(input("Enter time of last periapsis (in epoch timestamp) = "))
Pds=float(input("Enter period of orbits (in years) = "))
argperi_deg=float(input("Enter Argument of periapses in orbital plane (in degrees) = "))
m1=float(input("Enter mass of primary star (in terms of solar mass) = "))
m2=float(input("Enter mass of secondary star (in terms of solar mass) = "))
alpha1=float(input("Enter value of alpha in the luminosity relation for the primary star = "))
alpha2=float(input("Enter value of alpha in the luminosity relation for the secondary star = "))
#**************************************************************** TIME STAMP CONVERSION **********************************************

date_time1=datetime.utcfromtimestamp(t) 
time1=date_time1.strftime("%m/%d/%Y, %H:%M:%S")
#print(time1)
date_time2=datetime.utcfromtimestamp(T) 
time2=date_time2.strftime("%m/%d/%Y, %H:%M:%S")
#print(time2)

Time_period_of_obs=T-t
date_time_final=datetime.utcfromtimestamp(Time_period_of_obs) 
time_final=date_time_final.strftime("%X")
print('Time Period of Observation', time_final)

#**************************************************************** FUNCTION CALLS **********************************************

statement=bool(sys_dys>=0.01 and sys_dys<=50 and ang_arcsec>=0.01 and ang_arcsec<=2.0 and inc_deg>=0.001 and inc_deg<=179.999 and longAN_deg>=0.001 and longAN_deg<=179.999 and e>=0.001 and e<=0.95 and m1>=0.1 and m1<=10 and m2>=0.1 and m2<=m1 and Pds>=1 and Pds<=100)
if statement==True:
        ang_arcsec=ang_arcsec*pi/100/3600
        inc_deg=radians(inc_deg)
        longAN_deg=radians(longAN_deg)
        argperi_deg=radians(argperi_deg)
        a1,a2,v1r,v2r,L1,L2,pos_angle,TA,e_latest_last,M_deg,n=orbit(t, sys_dys, ang_arcsec, inc_deg, longAN_deg, e, T, Pds, argperi_deg, m1, m2, alpha1,alpha2)
        print('/n *********************** OUTPUTS *********************** /n')
        print('Semi-Major axis of Primary Star : ',a1.real,' (parsec)')
        print('Semi-Major axis  of Secondary Star : ',a2.real,' (parsec)')
        print('Peak Radial Velocity of Primary Star : ',v1r.real,' (km/s)')
        print('Peak Radial Velocity of Secondary Star : ',v2r.real,' (km/s)')
        print('Luminosity of Primary Star : ',L1,' (in terms of Solar Luminosity)')
        print('Luminosity of Secondary Star : ',L2,' (in terms of Solar Luminosity)')
        print('Position Angle of the Binary Star System in Image-plane : ',pos_angle,' (degrees)')
        print('True Anomaly of the Binary Star System : ',TA.real,' (degrees)')
        print('Eccentric Anomaly of the Binary Star System : ',e_latest_last.real,' (degrees)')
        print('Mean Anomaly of the Binary Star System : ',M_deg,' (degrees)')
        print('Mean Motion of the Binary Star System : ',n,' (radians per year)')
        
        
        V1,V2=radial_velocity(Pds, argperi_deg, m1, m2)
        pl.plot(arange(0,Pds-1),V1)
        pl.plot(arange(0,Pds-1),V2)
        pl.grid()
        pl.xlabel('Time (yrs)')
        pl.ylabel('Radial Velocity () ')
        pl.legend(['Primary Star','Secondary Star'],loc='upper right')
        pl.title('Radial Velocity Curve for Binary Star')
        pl.show()
        
        print('Do you have the peak wavelength of the star? ')
        num=int(input('Enter 1 if YES. Enter 2 if NO.'))
        
        if num==1:
            star_call=int(input('Enter 1 to display Star Type of Primary Star and Enter 2 to display Star Type of Secondary Star'))
            if star_call==1:
                D,Temp=StarType(L1)
                pl.plot(D,Temp)
                pl.grid()
                pl.title('Primary Star')
                pl.xlabel('Distance from Star (AU)')
                pl.ylabel('Temperature of the Star (Solar Temperature)')
                pl.legend('Distance vs Temperature Curve')
                pl.show()
            elif star_call==2:
                D,Temp=StarType(L2)
                pl.plot(D,Temp)
                pl.grid()
                pl.title('Secondary Star')
                pl.xlabel('Distance from Star (AU)')
                pl.ylabel('Temperature of the Star (Solar Temperature)')
                pl.legend('Distance vs Temperature Curve')
                pl.show()
            else:
                print('Wrong value entered. Try again later.')
        elif num==2:
            image_file=str(input('Enter the name of image file (with .jpg)'))
            ImgPlot(image_file)
        else:
            print("Error Encountered!")
        
        
        Simulate(a1.real,a2.real,e,inc_deg,100)
            
    
else:
    print("Inputs entered were not in the right range ")
    






