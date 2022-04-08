from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
from scipy import interpolate
from osgeo import gdal
import matplotlib as mpl
import numpy
import shapefile
import os
from datetime import datetime, date, time,timedelta

fontsizeV = 16


mpl.rcParams['font.size'] = 18.
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['axes.labelsize'] = 18.
mpl.rcParams['xtick.labelsize'] = fontsizeV
mpl.rcParams['ytick.labelsize'] = fontsizeV


cmptmp = plt.get_cmap("jet") 

# configurations
xmin = 115.3
dx = 0.02
Nx = 100
xmax = xmin + (Nx-1)*dx

ymin = 39.4
dy = 0.02
Ny = 80
ymax = ymin + (Ny-1)*dy

Nt = 10000
Nz = 1

length = 4
filename = "/data/Polyphemus/CASES/polyphemus/genDispersion/data/ground/LUC-glcf.bin"
names = ["Water", "Forest", "Grassland", "Cropland", "Urban"]




# configurations
endTime =  datetime(2017,6,01,20,0)
beginTime =  datetime(2017,01,01,8,0)
tmp = endTime-beginTime
tmp = tmp.total_seconds()

t = int(tmp/3600.0)
z = 1


datestr= endTime.strftime('Local time: %H:%M, %b %d, %Y')

tnum = Nx*Ny*Nz*(t-1) + Nx*Ny*(z-1)
tnum = tnum*length

# draw figure
fig = plt.figure(figsize=(12, 8))
	
# read results from binary file

f = open(filename, "rb")

data = numpy.fromfile(f, dtype=numpy.float32, count=Nx*Ny*14)

data = numpy.reshape(data, (14, Ny,Nx))

x = numpy.linspace(xmin, xmax, data.shape[2])
y = numpy.linspace(ymin, ymax, data.shape[1])
xnew = numpy.linspace(xmin, xmax, data.shape[2]*5)
ynew = numpy.linspace(ymin, ymax, data.shape[1]*5)
xx, yy = numpy.meshgrid(xnew, ynew)

ids = [[0, 0],[1, 9], [10, 10], [11, 11], [13, 13]]
for i in range(0,5,1):
	m = Basemap(llcrnrlon=xmin,
    	llcrnrlat=ymin,
   	urcrnrlon=xmax,
   	urcrnrlat=ymax,
   	resolution = None, 
   	projection = 'cyl',
    	suppress_ticks = False)
	sf = m.readshapefile("./BeijingForFigure", 'beijing',linewidth=1.5)

	for info in m.beijing_info:
		name = info['NAME99']
		name = name.split()
		plt.text(info['CENTROID_X']-0.08,info['CENTROID_Y'],name[0], fontname='Arial', fontsize=fontsizeV)
	data2 = 0
	for j in range(ids[i][0],ids[i][1]+1):
		ft = interp2d(x, y, data[j,:,:], kind='linear')
		data2 = data2 + ft(xnew, ynew)
		
	m.pcolormesh(xx, yy, data2,  cmap=cmptmp)
	plt.xlabel('Longitude (degree)')
	plt.ylabel('Latitude (degree)')
	plt.title(names[i])	
	cbar = plt.colorbar(extend='both')

	cbar.ax.set_title('Fraction in grid',fontname='Arial', fontsize=fontsizeV) 
	plt.savefig(names[i]+".png")
	plt.clf()
		















#f.close()
# read binary file

#plt.show()

