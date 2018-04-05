import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.constants import G
from matplotlib import animation


#code that reads in the data from the files
data1 = ascii.read('table3.dat', readme = 'ReadMe')
data2 = ascii.read('table5.dat', readme = 'ReadMe')

#assigning the values of radial velocity to these variables below 
#as well as time, one for data set 1
#and another one for data set 2

#data set 1
rad_vel1 = data1['RV1']
rad_vel2 = data1['RV2']
time = data1['HJD']

#data set 2
rad_vel11 = data2['RV1']
rad_vel22 = data2['RV2']
time1 = data2['HJD']

#code that will index our radial velocity so that we can zoom in on a region of interest
#in this case this will find for data set 1 the radial velocity of RV1 from when
#time is less than 54500

#initializer
i = 0
#list that will hold the index
time_segment = []

#while loop that loops through the condition required in this case that time[i] < 54500
while i < len(time):
	if time[i] < 54500:
		time_segment.append(i)
	i+=1
#assigning the value of index to the length of the list this will be used below to 
#only get the time and radial velocity values that we want
index = len(time_segment)


#code that cuts away unwanted values and only focuses on the region of interest 
rv1_segment = rad_vel1[:index]
time_adjusted = time[:index]


#this code is an exact copy of the above just that now we are looking at a different range
#and applying it to a different variable rad_vel2

#initializer
i = 0
#list that will hold the index
time_segment1 = []

#while loop that loops through the condition required in this case that time[i] < 54500
while i < len(time):
	if time[i] > 56000 and time[i] < 57000:
		time_segment1.append(i)
	i+=1
start = time_segment1[0]
end = time_segment1[-1]



rv2_segment = rad_vel2[start:end]
rv1_segment2 = rad_vel1[start:end]
time_seg = time[start:end]


#code that plots the zoomed in region for radial velocity 1
plt.title('Zoomed in ragion for RV1')
plt.scatter(time_adjusted, rv1_segment, c = 'blue')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()

plt.figure()

plt.title('Zoomed in ragion for RV1')
plt.scatter(time_seg, rv1_segment2, c = 'blue')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()


plt.figure()

plt.title('Zoomed in region for RV2')
plt.scatter(time_seg, rv2_segment, c = 'red')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()

#from this plot and other zoomed in regions we can calculate the orbital period in days by #looking for one period of a sine wave 

plt.figure()

#this code block plots the entirety of the data from data set 1
plt.title('Plot of the entire data set for data1')
plt.plot(time, rad_vel1,'.', label = 'Radial Velocity 1', color = 'blue')
plt.plot(time, rad_vel2,'.', label = 'Radial Velocity 2', color = 'red')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.legend()
plt.show()

plt.figure()

#this code block plots the entirety of the data from data set 2
plt.scatter(time1, rad_vel11, c = 'blue')
plt.scatter(time1, rad_vel22, c = 'red')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()

#making the mass function calculation (f)
def mass_func(period, RV, e):
	prefactor = (period * RV**3)/(np.pi * 2 * G)
	f = prefactor * (1 - e**2)**(3/2)
	return f


#####
#code block for animation of the binary system 
#####

#fig will be the variable that will hold the stuff inside the figure
fig = plt.figure()

#this code just sets uo our axes width inside the figure.
ax = plt.axes(xlim(0,10), ylim =(-5,5))

#This is the code for the plot element 
line, = ax.plot([], [], lw = 2)

###
#we need a way to calculate the position of the two stars relative to one another 
#From this we can use the distance calculated above from the average separation about 
#their center. So we need to be able to solve differential equations to model
#the change in position of the stars as a function of time. 
#Thus we will need to solve newtons laws as applied to planetary motion.
