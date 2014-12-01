from struct import unpack
import numpy as np

def read_long(filename):
	
	data = np.genfromtxt(filename,comments='#')
	x = data[:,3]
	y = data[:,4]	
	z = data[:,5]

	nscat = len(x)
	
	return x,y,z,nscat
	
