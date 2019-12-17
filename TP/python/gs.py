# Copyright 2009,  Lambert, Paletou, Josselin, Glorian
# This file is part of RTtools.
# RTtools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# RTtools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with RTtools.  If not, see <http://www.gnu.org/licenses/>.


import RT
#
#--- This is a custom modules! ask RT.py to fpaletou@irap.omp.eu in case

#--- warning when s=sold, call to s(new)=RT.formalGS(sold), (sold-snew)
from copy import deepcopy

from pylab import *
from numpy import *

#--- Generates semi-infite slab grid
dtau1=1.e-2
taumax=1.e4
npdec=5
#--- Collisional destruction probability, or (1-albedo)
eps=1.e-4
#--- Number of iterations for the iterative process
niter=150
#--- the above will be an INPUT

#--- The logarithmic opacity/spatial grid ---
fac = 10.**(1./float(npdec))
dtau=dtau1
t=0.
k=1
while (t < taumax):
    k=k+1
    t=t+dtau
    dtau=fac*dtau
npts=k
z=zeros((npts))
z[0]=0.
z[1]=dtau1
dtau=dtau1*fac
for k in arange(2,npts):
    z[k]=z[k-1]+dtau
    dtau=fac*dtau
print z

#--- Initial source function (Planck) ---
s=ones((npts))
b=ones((npts))

#--- Absorption coefficient set to unity throughout the slab ---
chi=ones((npts))

#--- Upper boundary condition (in general MU and FREQ dependent)
bcu=0.

#--- Lower boundary condition (in general MU and FREQ dependent)
bcd=1.

#--- Zero boundary conditions (for Lstar computation)
zerobc=0.

#--- Two-stream approximation/Eddington's approximation
mu=1./sqrt(3.)
wtdir=0.5
#--- In general, use your preferred quadrature (e.g., Gauss-Legendre)

#--- Analytical reference solution of Eddington
sedd=1. - (1.-sqrt(eps))*( exp(-sqrt(3.*eps)*z) )

#--------------------------------
#--- Compute diagonal operator  |
#--------------------------------
lstar=RT.formal(z,npts,chi,zerobc,zerobc,s,mu,wtdir,True)

shist=zeros((npts,niter))
relerr=zeros((niter))
dsk=zeros((npts))
te=zeros((niter))

#--- Relaxation parameter SOR: 1 < omega < 2, GS: omega=1.
omega=1.5

#-----------------------
#--- GS/SOR iterations |
#-----------------------
for i in arange(niter):
    sold=deepcopy(s)
    RT.formalGS(z,npts,chi,bcu,bcd,s,mu,wtdir,False,b,eps,lstar,omega)
    relerr[i]=max(abs(s-sold)/s)
    te[i]=max(abs(s-sedd)/sedd)
    shist[:,i]=s[:]
    print 'Iter: ',i, 'Max.rel.err.S: ',relerr[i], 'True.err: ', te[i]

print(s)
        
#--- SOR may generate negative values of S at places: set @ 10^(1/2) for plot
for nit in arange(niter):
    for k in arange(npts):
        if (shist[k,nit] < 0.):
#            shist[k,nit]=10.**(0.5)
            shist[k,nit]=-shist[k,nit]

figure(1)
plot(log10(z[1:npts]),log10(shist[1:npts,:]))
plot(log10(z[1:npts]),log10(sedd[1:npts]),'o')
xlabel(r"log($\tau$)",fontsize=18)
ylabel('log(S/B)',fontsize=18)

figure(2)
semilogy(relerr, '--k', linewidth=2)
semilogy(te,'k',linewidth=2)
xlabel('iteration number',fontsize=18)
ylabel('Max. relative correction and true error',fontsize=18)

show()

