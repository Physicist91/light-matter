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

from pylab import *
from numpy import *

def formal(z,nd,chi,bcu,bcd,s,csz,wtdir,evaldiag):
#
#-------------------------------------------------#
#     Short-characteristics formal solver         #
# 1D plane-parallel radiation transfer equation   #
#-------------------------------------------------#
#
#  written by fpaletou@ast.obs-mip.fr (July 2009)
#         now fpaletou@irap.omp.eu (made public Oct 2015)
#
#--- please cite Lambert, Paletou, Josselin & Glorian, Eur. J. Phys.
#--- (arXiv:1509.01158) for any further use!
#
#--- references: http://cdsads.u-strasbg.fr/abs/1987JQSRT..38..325O
#                http://adsabs.harvard.edu/abs/1994A%26A...285..675A
#                http://cdsads.u-strasbg.fr/abs/1995ApJ...455..646T
#                http://adsabs.harvard.edu/abs/2003A%26A...411..221C
#                http://cdsads.u-strasbg.fr/abs/2007JQSRT.103...57P
#    (in french) http://tel.archives-ouvertes.fr/tel-00332781/
#
#--- aknowledgements: to my 'sensei' Larry H. Auer, and my friends
#                     Loic Chevallier and Ludovick Leger                
#
#--- Inputs:  z......... spatial grid (nd)
#             nd........ # points in z
#             chi....... absorption coefficient (nd)
#             bcu....... lower boundary external radiation
#             bcd....... upper boundary external radiation
#             s......... source function
#             csz....... direction cosine of ray
#             wtdir..... angular quadrature weight
#             evaldiag.. BOOLEAN for Lstar OR Jbar computation
#
#--- Output:  jbar...... either Jbar OR Lstar (diagonal operator)
#
    jbar=zeros((nd))

    for imu in range(-1,3,2):
#
#--- starts boundary conditions stuff for mu>0 or mu<0
#
        if (imu < 0):
            k0=1
            k1=nd
            kdel=1
            if (evaldiag):
                xint=0.
            else:
                xint=bcu
#
#--- CAUTION: consider that, in general, xint(freq,dir), so that
#---   instead of this simple command-line, a double integration
#---   over 'freq' (inner) and 'dir' (outer) should be used
#
            jbar[k0-1]=jbar[k0-1] + wtdir*xint
#
#----------------------------------------------------|
#
        else:
            k0=nd
            k1=1
            kdel=-1
            if (evaldiag):
                xint=0.
            else:
                xint=bcd
            jbar[k0-1]=jbar[k0-1] + wtdir*xint
#
#--- ends boundary conditions contribution to Jbar
#
        for k in range(k0+kdel,k1+kdel,kdel):
            mu=imu*csz
            ku=k-kdel
            du=(z[ku-1]-z[k-1])/mu
            chiu=chi[ku-1]
            chi0=chi[k-1]
            su=s[ku-1]
            s0=s[k-1]
#
            if (k==k1):
#
#--- forces linear interpolation at boundary (d: undefined there)
#
                dd=du
                sd=2.*s0-su
                chid=chiu
            else:
                kd=k+kdel
                dd=(z[k-1]-z[kd-1])/mu
                sd=s[kd-1]
                chid=chi[kd-1]
#
#--- CAUTION: consider that, in general one may expect, s(freq) -
#---   eventually s(dir,freq) too - and chi(freq,k)
#
            dtu=0.5*(chiu+chi0)*du
            dtd=0.5*(chid+chi0)*dd
            exu=exp(-dtu)
#
            if (dtu <= 0.01):
#
#--- this may definitely be useful, believe me...
#
                w0=dtu*(1.-dtu/2.+dtu**2/6.-dtu**3/24.+dtu**4/120. \
                            -dtu**5/720.+dtu**6/5040.-dtu**7/40320. \
                            +dtu**8/362880.)
                w1=dtu**2*(0.5-dtu/3.+dtu**2/8.-dtu**3/30.+dtu**4/144. \
                               -dtu**5/840.+dtu**6/5760.-dtu**7/45360. \
                               +dtu**8/403200.)
                w2=dtu**3*(1./3.-dtu/4.+dtu**2/10.-dtu**3/36. \
                               +dtu**4/168.-dtu**5/960.+dtu**6/6480. \
                               -dtu**7/50400.+dtu**8/443520.)
            else:
                w0= 1.-exu
                w1= w0-dtu*exu
                w2= 2.*w1-dtu*dtu*exu
#
            psi0=w0+ (w1*(dtu/dtd-dtd/dtu)-w2*(1./dtd+1./dtu))/(dtu+dtd)
            psiu=(w2/dtu + w1*dtd/dtu)/(dtu+dtd)
            psid=(w2/dtd - w1*dtu/dtd)/(dtu+dtd)
#
            if (evaldiag):
                xint=psi0
            else:
                xint=xint*exu + psiu*su + psi0*s0 + psid*sd
#
#--- caution: consider that, in general, xint(freq,dir), so that
#---   instead of this simple command-line, a double integration
#---   over 'freq' (inner) and 'dir' (outer) should be used
#
            jbar[k-1]=jbar[k-1] + wtdir*xint
#
    return jbar

def formalGS(z,nd,chi,bcu,bcd,s,csz,wtdir,evaldiag,b,eps,lstar,omega):
#
#-------------------------------------------------#
#     Short-characteristics formal solver         #
#     --- for GS/SOR iterative scheme ---         #
# 1D plane-parallel radiation transfer equation   #
#-------------------------------------------------#
#
#  written by fpaletou@ast.obs-mip.fr (July 2009)
#         now fpaletou@irap.omp.eu (made public Oct 2015)
#
#--- please cite Lambert, Paletou, Josselin & Glorian, Eur. J. Phys.
#--- (arXiv:1509.01158) for any further use!
#
#--- references: http://cdsads.u-strasbg.fr/abs/1987JQSRT..38..325O
#                http://cdsads.u-strasbg.fr/abs/1995ApJ...455..646T
#                http://cdsads.u-strasbg.fr/abs/2007JQSRT.103...57P
#    (in french) http://tel.archives-ouvertes.fr/tel-00332781/
#
#--- aknowledgements: to my 'sensei' Larry H. Auer, and my friends
#                     Loic Chevallier and Ludovick Leger                
#
#--- Inputs:  z......... spatial grid (nd)
#             nd........ # points in z
#             chi....... absorption coefficient (nd)
#             bcu....... lower boundary external radiation
#             bcd....... upper boundary external radiation
#             s......... source function
#             csz....... direction cosine of ray
#             wtdir..... angular quadrature weight
#             evaldiag.. BOOLEAN for Lstar OR Jbar computation
#             eps....... collisional destruction parameter [GS]
#             lstar..... diagonal operator [GS]
#             omega..... relaxation parameter [GS]
#
#--- Output:  S......... S(new) of GS/SOR
#
    jbar=zeros((nd))
 #---GS specific
    psidd=zeros((nd))

    for imu in range(-1,3,2):
#
#--- starts boundary conditions stuff for mu>0 or mu<0
#
        if (imu < 0):
            k0=1
            k1=nd
            kdel=1
            if (evaldiag):
                xint=0.
            else:
                xint=bcu
#
#--- CAUTION: consider that, in general, xint(freq,dir), so that
#---   instead of this simple command-line, a double integration
#---   over 'freq' (inner) and 'dir' (outer) should be used
#
            jbar[k0-1]=jbar[k0-1] + wtdir*xint
#
#----------------------------------------------------|
#
        else:
            k0=nd
            k1=1
            kdel=-1
            if (evaldiag):
                xint=0.
            else:
                xint=bcd
            jbar[k0-1]=jbar[k0-1] + wtdir*xint
#
#--- ends boundary conditions contribution to Jbar


#---GS specific: update S at boundaries, before advancing to next k
            dsk=((1.-eps)*jbar[k0-1]+eps*b[k0-1]-s[k0-1]) \
                /(1.-(1.-eps)*lstar[k0-1])
            s[k0-1]=s[k0-1] + omega*dsk
#
#----------------------------------------------------|

        for k in range(k0+kdel,k1+kdel,kdel):
            mu=imu*csz
            ku=k-kdel
            du=(z[ku-1]-z[k-1])/mu
            chiu=chi[ku-1]
            chi0=chi[k-1]
            su=s[ku-1]
            s0=s[k-1]
#
            if (k==k1):
#
#--- forces linear interpolation at boundary (d: undefined there)
#
                dd=du
                sd=2.*s0-su
                chid=chiu
            else:
                kd=k+kdel
                dd=(z[k-1]-z[kd-1])/mu
                sd=s[kd-1]
                chid=chi[kd-1]
#
#--- CAUTION: consider that, in general one may expect, s(freq) -
#---   eventually s(dir,freq) too - and chi(freq,k)
#
            dtu=0.5*(chiu+chi0)*du
            dtd=0.5*(chid+chi0)*dd
            exu=exp(-dtu)
#
            if (dtu <= 0.01):
#
#--- this may definitely be useful, believe me...
#
                w0=dtu*(1.-dtu/2.+dtu**2/6.-dtu**3/24.+dtu**4/120. \
                            -dtu**5/720.+dtu**6/5040.-dtu**7/40320. \
                            +dtu**8/362880.)
                w1=dtu**2*(0.5-dtu/3.+dtu**2/8.-dtu**3/30.+dtu**4/144. \
                               -dtu**5/840.+dtu**6/5760.-dtu**7/45360. \
                               +dtu**8/403200.)
                w2=dtu**3*(1./3.-dtu/4.+dtu**2/10.-dtu**3/36. \
                               +dtu**4/168.-dtu**5/960.+dtu**6/6480. \
                               -dtu**7/50400.+dtu**8/443520.)
            else:
                w0= 1.-exu
                w1= w0-dtu*exu
                w2= 2.*w1-dtu*dtu*exu
#
            psi0=w0+ (w1*(dtu/dtd-dtd/dtu)-w2*(1./dtd+1./dtu))/(dtu+dtd)
            psiu=(w2/dtu + w1*dtd/dtu)/(dtu+dtd)
            psid=(w2/dtd - w1*dtu/dtd)/(dtu+dtd)
#
#---GS specific: store those Psi's reused for add. corrections
            if (imu < 0):
                psidd[k-1]=psid
            else:
                psi0u=psi0
#---------------------------------|
            if (evaldiag):
                xint=psi0
            else:
                xint=xint*exu + psiu*su + psi0*s0 + psid*sd
#
#--- caution: consider that, in general, xint(freq,dir), so that
#---   instead of this simple command-line, a double integration
#---   over 'freq' (inner) and 'dir' (outer) should be used
#

#---GS specific: upward Delta Jbar correction
            if ((imu<0) or (evaldiag)):
                jbar[k-1]=jbar[k-1] + wtdir*xint
            else:
#---GS specific: Eq. (39) of Trujillo Bueno & Fabiani Bendicho (1995)
                jbar[k-1]=jbar[k-1] + wtdir*xint \
                    + omega*dsk*wtdir*psidd[k-1]
#
#---GS specific: AFTER inner freq and dir loops, in general----|
#
            if ((imu>0) and (not evaldiag)):
                dsk=((1.-eps)*jbar[k-1]+eps*b[k-1]-s[k-1]) \
                    /(1.-(1.-eps)*lstar[k-1])
                s[k-1]=s[k-1] + omega*dsk
#---GS specific: Eq. (40) of Trujillo Bueno & Fabiani Bendicho (1995)
                xint=xint + omega*dsk*psi0u
                
#---GS specific: now return S, not jbar - S(new) now exists!
    return s
