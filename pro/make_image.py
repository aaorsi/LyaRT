from read_LyaRT import *
from read_LyaRT_long import *
from matplotlib import pyplot as pl
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
import os


filename='/home/aaorsi/work/LyaRT/Output/Bicone_test/test_bicone'
wl, wl0,x,y,z,theta,phi,Np,Nout,nscat, interact=read_wave(filename+'.data')

nphotons = 300
nscatmax = 2e3

ylim = 2e20
zlim = 2e20
xlim = 2e20


x_s = np.zeros([nscatmax,nphotons])
y_s = np.zeros([nscatmax,nphotons])
z_s = np.zeros([nscatmax,nphotons])


i = 0
for j in np.arange(nphotons):
	print j
	FileLong = filename+'_long.data'+str(j)
	if os.path.isfile(FileLong):
		x_, y_, z_, iscat = read_long(FileLong)
		x_s[:iscat,i] = x_/zlim
		y_s[:iscat,i] = y_/zlim
		z_s[:iscat,i] = z_/zlim
		i += 1
	


app_allowed = 25.0   # degrees

app = app_allowed * np.pi/180.0

view_angle = np.pi/2.

r_inner = 5.0e19
apperture = np.pi/4.

ip = np.where((theta >= view_angle - app) &
			  (theta <= view_angle + app) ) 	

nx = 100.0
ny = 100.0
nz = 100.0


#ymay = np.may(y[ip])
#ymin = np.min(y[ip])
#zmay = np.may(z[ip])
#zmin = np.min(z[ip])

ymax = ylim
ymin = -ylim
zmax = zlim
zmin = -zlim
xmax = xlim
xmin = -xlim

#y[ip] = 2*((y[ip]-ymin)/(ymay-ymin))-1
#z[ip] = 2*((z[ip]-zmin)/(zmay-zmin))-1

ybin = (ymax-ymin)/(ny-1)
zbin = (zmax-zmin)/(nz-1)
xbin = (xmax-xmin)/(nx-1)


yz_image = np.zeros([ny,nz])


for iy in np.arange(ny):
	_y = iy * ybin + ymin
	for iz in np.arange(nz):
		_z = iz *zbin + zmin
		sel = np.where((y[ip] >= _y-ybin/2.) &
					   (y[ip] <= _y+ybin/2.) &
					   (z[ip] >= _z-zbin/2.) &
					   (z[ip] <= _z+zbin/2.))

#		print 'iy=',_y,'iz=',_z,'len=',len(sel[0])

		yz_image[iy][iz] = len(sel[0])



a1 =  pl.figure()
ax = a1.add_axes([0,0,1,1],projection='3d')
ax.axis('off')
#ax = pl.axes(xlim=(-1,1),ylim=(-1,1))
#line, = ax.plot([],[],'-y',lw=1)

colors = pl.cm.jet(np.linspace(0, 1, nphotons))
lines = sum([ax.plot([],[],[],'-', c=c)
		for c in colors],[])
pts = sum([ax.plot([],[],[],'o',c=c,markersize=2.5)
		for c in colors],[])


#pl.imshow(yz_image,extent=[-1,1,-1,1],interpolation='nearest',origin='lower',alpha=0.70,cmap='gray')
#pl.xticks(np.arange(-1,1,5))
#pl.yticks(np.arange(-1,1,5))
#pl.xlabel('y')
#pl.ylabel('z')    
#pl.autoscale(False)
#cbar = pl.colorbar()


ax.set_xlim((-1,1))
ax.set_ylim((-1,1))
ax.set_zlim((-1,1))

ax.view_init(30,0)

theta = np.arange(-apperture,apperture,0.01)
rad = r_inner/zlim
rx = rad * np.cos(theta)
ry = rad * np.sin(theta)


#rad = theta * 0 + r_inner/xlim
#pl.plot(ry,rx,'--r')
#pl.plot(-ry,-rx,'--r')
#pl.plot(3*ry,3*rx,'-.r')
#pl.plot(-3*ry,-3*rx,'-.r')
#rr = np.arange(r_inner/zlim,3*r_inner/zlim,0.1)
#th0 = apperture

#rx = rr * np.cos(th0)
#ry = rr * np.sin(th0)
#pl.plot(ry,rx,'--r')
#pl.plot(-ry,-rx,'--r')

#th0 = -apperture

#rx = rr * np.cos(th0)
#ry = rr * np.sin(th0)
#pl.plot(ry,rx,'--r')
#pl.plot(-ry,-rx,'--r')

nframes = 600

def init():
	for line,pt in zip(lines, pts):
		line.set_data([],[])
		line.set_3d_properties([])

		pt.set_data([],[])
		pt.set_3d_properties([])
	return lines + pts


def animate(i):	
	global nphotons, x_s,y_s,z_s
	ip = 0
	for line, pt  in zip(lines,pts):
		insc = np.where(y_s[:,ip] != 0)
		nsc = len(insc[0])
		_i = nsc-1 if i >= nsc else i
#		if i == 0:
#			line.set_data([0,0], [0,0])
#			line.set_3d_properties([0,0])
#
#			pt.set_data([0],[0])
#			pt.set_3d_properties([0])
#		else:
		line.set_data(x_s[:_i,ip],y_s[:_i,ip])
		line.set_3d_properties(z_s[:_i,ip])
		pt.set_data(x_s[_i,ip],y_s[_i,ip])
		pt.set_3d_properties(z_s[_i,ip])			

		ip += 1
#		print i, ip	
			

	ax.view_init(30+0.01*i,0.3*i)
	a1.canvas.draw()
#	print y_s[_i,ip],z_s[_i,ip]
	return lines + pts

anim = animation.FuncAnimation(a1,animate,init_func = init,frames=nframes,blit=True)



#pl.plot(theta,-rad,'r',linewidth=3)
#pl.plot(theta,rad*10,'r',linewidth=3)
#pl.plot(theta,-rad*10,'r',linewidth=3)
#r0 = np.arange(r_inner/xlim,r_inner*10/xlim,0.1)
#th0 = r0*0 + apperture
#pl.plot(th0,r0,'r',linewidth=3)
#pl.plot(-th0,r0,'r',linewidth=3)
#pl.plot(th0,-r0,'r',linewidth=3)
#pl.plot(-th0,-r0,'r',linewidth=3)
#


anim.save('test.mp4',fps=25,bitrate=-1)

#pl.show()

#pl.savefig('cone_faceon.png',bbox_inches='tight',dpi=350)
pl.close(a1)
#
