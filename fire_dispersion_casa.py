# PROBLEM 2
#
# This program simulates a wildfire in a square area of forest.  Taking into 
# account diffusion, heat loss, wind, and combustion, use the Forward Euler
# Method to show how the temperature and wood density at a given location change
# with time.  Then set up an initial distribution for the wood density.  Pretend 
# that the square of forest has infinitesimally thin parallel lines of constant wood
# density.  These lines all have the same slope, so their y-intercepts can serve
# as indicators of their respective wood densities.  Please see the question
# introduction video for more details.

import math
import matplotlib.pyplot
import numpy 
from libtiff import TIFF

diffusion_coefficient = 5. # m2 / s
ambient_temperature = 310. # K
heat_loss_time_constant = 120. # s
velocity_y = 0.4 # m / s
velocity_x = -0.2 # m / s
ignition_temperature = 561. # K valor inicial: 561
burn_temperature = 1400. # K
burn_time_constant = 0.5 * 3600. # s
heating_value = (burn_temperature - ambient_temperature) / (
heat_loss_time_constant * 100.) * burn_time_constant # K / (kg / m2)
slope = 0.4 # dimensionless
intercept_1 = 100. # m
intercept_2 = 170. # m
wood_1 = 100. # kg / m2
wood_2 = 70. # kg / m2

diagonal_factor = .707

length = 650. # meters; domain extends from -length to +length
# A grid size of 50 x 50 ist much too small to see the correct result. For a better result, set the size to 200 x 200. That computation would, however, be far too long for the Web-based development environment. You may want to run it offline.
size = 120 # number of points per dimension of the grid
dx = 2. * length / size
# Pick a time step below the threshold of instability
h = 0.2 * dx ** 2 / diffusion_coefficient # s
end_time = 30. * 60. # s

def bezier(p, t):
	# calculo de los coeficientes polinomiales
	#~ cx = 3 * (p[1][0] - p[0][0])
	#~ bx = 3 * (p[2][0] - p[1][0])
	#~ ax = p[3][0] - p[0][0] - cx - bx
	
	#~ cy = 3 * (p[1][1] - p[0][1])
	#~ by = 3 * (p[2][1] - p[1][1])
	#~ ay = p[3][1] - p[0][1] - cy - by
	
	result = []
	
	result.append( ( 1 - t ) ** 3 * p[0][0] + 3 * ( 1 - t ) ** 2 * t * p[1][0] + 3 * ( 1 -t ) * t ** 2 * p[2][0] + t ** 3 * p[3][0] )
	result.append( ( 1 - t ) ** 3 * p[0][1] + 3 * ( 1 - t ) ** 2 * t * p[1][1] + 3 * ( 1 -t ) * t ** 2 * p[2][1] + t ** 3 * p[3][1] )
	
	return result

# Convert from integer grid positions to coordinates measured in meters
def grid2physical(i, j):
	return i * dx - length + 0.5 * dx, j * dx - length + 0.5 * dx

def wildfire():
	temperatures_old = ambient_temperature * numpy.ones([size, size]) # K
	data = numpy.zeros((size, size), dtype = 'uint16')
	wood_old = 1000 * numpy.ones([size, size]) # kg / m2
	
	for i in range(40):
		pto1 = [[82,118],[89,98],[97,85],[118,78]]
		p = bezier(pto1, float(i+1) / 40)
		temperatures_old[int(round(p[0]))][int(round(p[1]))] = burn_temperature
            
	#for i in range(size):
	#	for j in range(size-3, size ):
			#~ x, y = grid2physical(i, j)
	#		temperatures_old[j][i] = burn_temperature
			# Given:
			# wood_old[j][i] = 100.
			
			# Your code here
	
	temperatures_new = numpy.copy(temperatures_old) # K
	wood_new = numpy.copy(wood_old) # kg / m2
		
	num_steps = int(end_time / h)

	for step in range(num_steps):
		
		if step <  size / 2:
			pto1 = [[38,30],[38,30],[50,53],[50,53]]
			t = float(step) / (size / 2)
			p = bezier(pto1, t)
			temperatures_old[int(round(p[0]))][int(round(p[1]))] = burn_temperature
                        

		if step > size / 2 and step < size:
			pto4 = [[50,53],[50,53],[90,46],[90,46]]
			t = float(step - (size / 2)) / (size / 2)
			p2 = bezier(pto4, t)
			temperatures_old[int(round(p2[0]))][int(round(p2[1]))] = burn_temperature

		if step > size and step < 3 * size / 2:
			pto4 = [[90,46],[90,46],[88,15],[88,15]]
			t = float(step - size) / (size / 2)
			p2 = bezier(pto4, t)
			temperatures_old[int(round(p2[0]))][int(round(p2[1]))] = burn_temperature

		#~ if step > 3 * size and step < 4 * size:
			#~ pto4 = [[96,22],[103,23],[105,26],[107,32]]
			#~ t = float(step - 3 * size) / size
			#~ p2 = bezier(pto4, t)
			#~ temperatures_old[int(round(p2[0]))][int(round(p2[1]))] = burn_temperature

		for j in range(1, size - 1):
			for i in range(1, size - 1):
				burn_rate = 0
				if temperatures_old[j][i] > ignition_temperature:
					burn_rate = wood_old[j][i] / burn_time_constant
				wood_new[j][i] = wood_old[j][i] - h * burn_rate
				temp = temperatures_old[j][i]
				temperatures_new[j][i] = temp + h * (diffusion_coefficient / dx ** 2 * (temperatures_old[j][i+1] + temperatures_old[j][i-1] + temperatures_old[j+1][i] + temperatures_old[j-1][i] + diagonal_factor * (temperatures_old[j+1][i+1] + temperatures_old[j+1][i-1] + temperatures_old[j-1][i+1] + temperatures_old[j-1][i-1])- (4 + 4 * diagonal_factor) * temp) - 0.5 / dx * (velocity_x * (temperatures_old[j+1][i] - temperatures_old[j-1][i]) + velocity_y * (temperatures_old[j][i+1] - temperatures_old[j][i-1]))) - (temp - ambient_temperature) / heat_loss_time_constant + heating_value * burn_rate
				
				if data[j][i] == 0:
					if temperatures_new[j][i] > ignition_temperature:
						data[j][i] = step * 65530 / (num_steps)
		temperatures_old, temperatures_new = temperatures_new, temperatures_old
		wood_old, wood_new = wood_new, wood_old

	img = TIFF.open('p3.tiff', mode='w')
	img.write_image(data)
	img.close()
	return temperatures_old


temperatures_old = wildfire()

def fire_plot():
	dimensions = [-length, length, -length, length]
	 
	axes = matplotlib.pyplot.subplot(121)
	matplotlib.pyplot.imshow(temperatures_old, interpolation = 'nearest', cmap = matplotlib.cm.hot, origin = 'lower', extent = dimensions)
	matplotlib.pyplot.colorbar()
	axes.set_title('Temperature in K')
	axes.set_xlabel('x in m')
	axes.set_ylabel('y in m')

#~ axes = matplotlib.pyplot.subplot(122)
#~ matplotlib.pyplot.imshow(wood_old, cmap = matplotlib.cm.winter, origin = 'lower', extent = dimensions)
#~ matplotlib.pyplot.colorbar()
#~ axes.set_title('Density of wood in kg/m$^2$')
#~ axes.set_xlabel('x in m')
#~ axes.set_ylabel('y in m')

fire_plot()    
