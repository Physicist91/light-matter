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
#-- This is a custom module! ask RT.py to fpaletou@irap.omp.eu in case

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
s=ones(npts)
b=ones(npts)

#--- Absorption coefficient set to unity throughout the slab ---
chi=ones(npts)

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

shist=zeros((npts,niter), 'Float64')
relerr=zeros((niter), 'Float64')
dsk=zeros((npts), 'Float64')
te=zeros((niter), 'Float64')

#------------------------------------
#--- Picard (aka. Lambda-iteration) |
#------------------------------------
for i in arange(niter):
    jbar=RT.formal(z,npts,chi,bcu,bcd,s,mu,wtdir,False)
    shist[:,i]=(1.-eps)*jbar + eps*b
    relerr[i]=max(abs(s-shist[:,i])/s)
    te[i]=max(abs(sedd-shist[:,i])/sedd)
    print 'Iter: ',i, 'Max.rel.err.S: ',relerr[i], 'True.err: ', te[i]
    s=shist[:,i]

print(s)
    
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

