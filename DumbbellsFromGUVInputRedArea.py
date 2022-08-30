########################################################################
# Interplay of membrane curvature and actin assembly                   #
# Modeling the shape and size of GUV dumbbells                         #
# Felix Frey, TU Delft, 2022                                           #   
########################################################################


#import used packages

import numpy as np
import pandas as pd
from scipy.optimize import minimize


#function definitions

#solve the volume equation for the opening angle of the bright spherical cap
#of the dumbbell theta_1 by minimization
#input: GUV volume V, surface area of the internal vesilce A_v, 
#total surface area of the GUV A_tot
def Theta_1_minimization(x,V,A_v,A_total):
    R_1=Radius_1(x,A_v)
    theta_2=Theta_2(x,A_v,A_total)
    R_2=Radius_2(theta_2,A_v,A_total)
    return (V-np.pi/3*R_1**3*(2+np.cos(x))*(1-np.cos(x))**2-np.pi/3*R_2**3*(2+np.cos(theta_2))*(1-np.cos(theta_2))**2)**2    
   
#opening anlge of the dark spherical cap theta_2
def Theta_2(theta_1,A_v,A_total):
    argument=(1+np.cos(theta_1))*A_v/(A_total-A_v)-1
    return np.arccos(argument)   

#radius of the bright spherical cap
def Radius_1(theta_1,A_v):
    return np.sqrt(A_v/(2*np.pi*(1-np.cos(theta_1))))

#radius of the dark spherical cap
def Radius_2(theta_2,A_v,A_total):
    return np.sqrt((A_total-A_v)/(2*np.pi*(1-np.cos(theta_2))))

#fit linear function
def func(x, m, b):
    return m*x+b

#write data to file
def WriteSimulationData(diameter_1_list,diameter_2_list,average_diameter,neck_diameter):
    np.savetxt("Simulated_Data.csv", np.c_[diameter_1_list,diameter_2_list,average_diameter,neck_diameter],header="bright half diameter (µm), dark half diameter (µm), average diameter (µm), neck diameter (µm)", delimiter=',')


#read in data from the measured GUV diameter distribution


data_file = pd.read_csv('GUVs_sizes.csv',sep=',')
processed = ["number","diameter (um)"]
df = data_file[processed]
df=df.to_numpy()
number=df[:,0]
diameter=df[:,1]


#parameter values for sampling


#number of simulated GUVs    
guv_number=10000

#mean and standard deviation of the Gaussian distribution for the reduced area
mean_reduced_area_GUV=1
std_reduced_area_GUV=0.03

#max GUV diameters
max_GUV_diameter=50

#counting variable
i=0

# data containers
theta_1_list=[]
theta_2_list=[]
diameter_1_list=[]
diameter_2_list=[]
average_diameter=[]
neck_diameter=[]


########################################################################
#main function
########################################################################

while(True):
  
    #check for the maximal number of GUVs to end the sampling
    if i>=guv_number:
        break
    
    #sample a reduced area and check if it is larger than 1.
    nu_large=np.random.normal(mean_reduced_area_GUV,std_reduced_area_GUV)
    nu_small=np.random.normal(mean_reduced_area_GUV,std_reduced_area_GUV)
    if  nu_large>=1 and nu_small>=1:
        
        #sample two radii
        radius_large = np.random.choice(diameter)/2
        radius_small=np.random.choice(diameter)/2
    
        # ensure to save the larger radius as radius and the smaller radius as
        if radius_large<radius_small:
            temp=radius_small
            radius_small=radius_large
            radius_large=temp
    
        #calculate the GUV volume, the GUV total surface area A_total and
        #surface area of the internal vesilce A_v
        V_large=4/3*np.pi*radius_large**3
        A_total=4*np.pi*radius_large**2*nu_large
        A_v=4*np.pi*radius_small**2*nu_small

        #solve the volume equation for the opening angle of the bright spherical cap theta_1
        #i.e. V(theta_1)=V_large, to find theta_1    
        #as an initial value we use pi/4, which is in the middle of the considered interval
        x_initial_0=np.pi/4
        root = minimize(Theta_1_minimization, x_initial_0,args=(V_large,A_v,A_total),method='Nelder-Mead')
          
        #save data
        theta_1=root.x[0]
        theta_2=Theta_2(theta_1,A_v,A_total)
        radius_1=Radius_1(theta_1,A_v)
        radius_2=Radius_2(theta_2,A_v,A_total)
                    
        #check if theta_1 and theta_2 are within the bounds (0,pi/2) and if the diameters are smaller than max_GUV_diameter
        if theta_2>0 and theta_2<np.pi/2 and theta_1>0 and theta_1<np.pi/2 and radius_1<max_GUV_diameter/2 and radius_2<max_GUV_diameter/2:
        
            #save data
            theta_1_list.append(theta_1)
            theta_2_list.append(theta_2)
            diameter_1_list.append(radius_1*2)
            diameter_2_list.append(radius_2*2)
        
            average_diameter.append((radius_1+radius_2))
            neck_diameter.append(2*np.sin(theta_1)*radius_1)
            i+=1

#write simulation data to file
WriteSimulationData(diameter_1_list,diameter_2_list,average_diameter,neck_diameter)







