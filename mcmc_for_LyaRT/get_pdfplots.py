# script to plot the MCMC steps, do contours and plot the PDFs

import pylab as pl
import numpy as np

datan = '2'
data =  np.loadtxt('spec_mcmc_new_'+datan+'.out') 
print len(data[:,0]), 'steps:' + np.str(len(data[:,0])/60.0) 
x, y = data[:,0], data[:,1]

nd = len(x) 

xmax, xmin = 20, 18
ymax, ymin = 100, 0

#xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()

nbins = 50.0
xbsize = (xmax - xmin)/(nbins -1)
ybsize = (ymax - ymin)/(nbins -1)

xarr = np.linspace(xmin,xmax,nbins)
yarr,step = np.linspace(ymin,ymax,nbins,retstep=True)

print step,ybsize

dgrid = np.zeros([nbins,nbins])

for j in range(nd):
	ix = int((x[j] - xmin) / xbsize)
	iy = int((y[j] - ymin) / ybsize)

	if ix < 0 or ix > nbins-1 or iy < 0 or iy > nbins-1:
		continue

	dgrid[iy][ix] += 1


#pl.contour(xarr,yarr,dgrid,10,origin='upper')

#pl.plot(data[:nd-60,0],data[:nd-60,1],'.') 
#pl.plot(x[nd-60:],y[nd-60:],'o') 
pl.xlabel(r'$\log(N_H)$',fontsize=20)
pl.ylabel(r'$V_{\rm exp}$',fontsize=20) 



pl.contour(xarr,yarr,dgrid,4,colors='gray',linewidths=2)

mesh = pl.pcolormesh(xarr,yarr,dgrid,cmap='YlGnBu',alpha=0.75)
pl.axis([xmin,xmax,ymin,ymax])


pl.savefig(datan + '_pdf.eps')
#pl.show()



