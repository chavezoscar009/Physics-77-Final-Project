import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.constants import G
from matplotlib import animation


#code that reads in the data from the files
data1 = ascii.read('table3.dat', readme = 'ReadMe')
data2 = ascii.read('table5.dat', readme = 'ReadMe')

#first step we need to find the period of the binary system 
#we do this by plotting the radial velocities and then we will zoom in on
#a region and see if we can deduce the period from it. 
#once we have the period all we need next will be to find out the mass and average
#separation 

#we can find ratio of mass by taking the ratio of radial velocities
 

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


zoomRV1 = []
i=0
while i < len(time):
	if time[i] > 56200 and time[i] < 56260:
		zoomRV1.append(i)
	i+=1
zoomstart = zoomRV1[0]
zoomend = zoomRV1[-1]

RV1_zoom = rad_vel1[zoomstart : zoomend]
timezoom = time[zoomstart : zoomend]


zoomRV2 = []
i=0
while i < len(time):
	if time[i] > 56200 and time[i] < 56300:
		zoomRV2.append(i)
	i+=1
zoomstart2 = zoomRV2[0]
zoomend2 = zoomRV2[-1]

RV2_zoom = rad_vel2[zoomstart2 : zoomend2]
tprime = time[zoomstart2 : zoomend2]

rv2_segment = rad_vel2[start:end]
rv1_segment2 = rad_vel1[start:end]
time_seg = time[start:end]


#code that plots the zoomed in region for radial velocity 1
plt.title('Zoomed in region for RV1')
plt.scatter(time_adjusted, rv1_segment, c = 'blue')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()

plt.figure()

plt.title('Zoomed in region for RV1')
plt.scatter(time_seg, rv1_segment2, c = 'blue')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()


plt.figure()

#code that plots the zoomed in region for radial velocity 1
plt.title('Zoomed in region for RV1')
plt.scatter(timezoom, RV1_zoom, c = 'blue')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity (km/s)')
plt.savefig('Rad_vel1.png')
plt.show()

plt.figure()

plt.title('Zoomed in region for RV2')
plt.scatter(time_seg, rv2_segment, c = 'red')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity')
plt.show()

plt.figure()

plt.title('Zoomed in region for RV2')
plt.scatter(tprime, RV2_zoom, c = 'red')
plt.xlabel('Time (days)')
plt.ylabel('Radial Velocity (km/s)')
plt.xlim(min(tprime), max(tprime))
plt.savefig('Rad_vel2.png')
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


###
#from the plots above we are able to infer that the orbital period for star 1 is of order #35 days +/- 2 days
#and for the second object/star, we can infer that the period is about 35 days as well.
###

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
#ax = plt.axes(plt.xlim(0,10), plt.ylim =(-5,5))

#This is the code for the plot element 
#line, = ax.plot([], [], lw = 2)

###
#we need a way to calculate the position of the two stars relative to one another 
#From this we can use the distance calculated above from the average separation about 
#their center. So we need to be able to solve differential equations to model
#the change in position of the stars as a function of time. 
#Thus we will need to solve newtons laws as applied to planetary motion.

###
#writing down the equation of planetary motion that we will be implementing
#I will be writing them as function that we can call and if we need to increment 
#these as a function of time then only time will vary
###


#####
#variables that we will use
#
#r12 will be the distance away between the two objects/stars
#x1 will represent the position in the x axis of object/star 1
#x2 will represent the position in the x axis of object/star 2
#y1 will represent the position in the y axis of object/star 1
#y1 will represent the position in the y axis of object/star 2
#v_x1 will represent the velocity in the x axis of object/star 1
#v_x2 will represent the velocity in the x axis of object/star 2
#v_y1 will represent the velocity in the y axis of object/star 1
#v_y2 will represent the velocity in the y axis of object/star 2
#m1 will represent the mass of star 1
#m2 will represent the mass of star 2


#for this will employ the center of mass reference frame which will give us our initial
#conditions necessary to at least set up the problem/orbit

#with COM we have:
#x1(0) = r0
#y1(0) = 0
#x2(0) = -(m1/m2)r0
#y2(0) = 0
#v_x1 = 0
#v_y1 = v0
#v_x2 = 0
#v_y2 = -(m1/m2)r0



#however since we have radial velocity data we can substitute these values in for v0
#this will be where the radial velocity is a maximum for one star

#where r0 is the elliptical path of the orbit given by:
#r0 = (1-e)a / (1+(m1/m2))
#and v0 is given by:
#(1/1+(m1/m2)) * sqrt((1+e)/(1-e)) *sqrt(Gm/a)

#here m is the total mass of the system m1+m2 and a is the semi-major axis
#and e is the eccentricity
#for circular orbits the eccentricity is 0




#function that calculates the separation between the two objects/stars
def separation(x1,x2,y1,y2):
	r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
	return r12

#function that calculates the position of x1 as a function of time
#for the functions below t is actually delta t, so a small time interval

def x1(x1, v_x1, t):
	new_x = x1 + v_x1*t
	return new_x

#function that calculates position of x2 as a function of time
def x2(x2, v_x2, t):
	new_x = x2 + v_x2*t
	return new_x

#function that calculates the position of y1 as a function of time
def y1(y1, v_y1, t):
	new_y = y1 + v_y1*t
	return new_y

#function that calculates the position of y2 as a function of time
def y2(y2, v_y2, t):
	new_y = y2 + v_y2*t
	return new_y

#below I will write the functions that will calculate the velocity components 
def vx1(x1,x2,r12,t,v_x1):
	a = (G * m2 * (x2-x1)) / r12**3
	new_vx1 = v_x1 + t*a
	return new_vx1

def vx2(x1,x2,r12,t,v_x2):
	a = (G * m1 * (x1-x2)) / r12**3
	new_vx2 = v_x2 + t*a
	return new_vx2


def vy1(y1,y2,r12,t,v_y1):
	a = (G * m2 * (y2-y1)) / r12**3
	new_vy1 = v_y1 + t*a
	return new_vy1

def vy2(y1,y2,r12,t,v_y2):
	a = (G * m1 * (y1-y2)) / r12**3
	new_vy2 = v_y2 + t*a
	return new_vy2




















