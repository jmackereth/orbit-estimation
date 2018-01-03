import os, os.path
import sys
import pickle
import numpy 
import numpy as np
from galpy.potential import MWPotential2014, rl
from galpy.potential import evaluatePotentials as evalPot
from galpy.orbit import Orbit
from galpy.actionAngle import estimateDeltaStaeckel
from galpy.actionAngle import actionAngleStaeckel, UnboundError
from galpy.actionAngle_src.actionAngleStaeckel import actionAngleStaeckelSingle, estimateDeltaStaeckel
from galpy.util import bovy_coords
from astropy import units


def orbit_at_E_L(logE, logL, x=4./5., pot=MWPotential2014):
    '''
    return an orbit instance at a given energy and angular momentum in a given
    potential
    '''
    Einf= evalPot(pot,10.**12.,0.)
    Rc= rl(pot,10**logL)
    Ec= evalPot(pot,Rc,0.)+0.5*(10**logL)**2./Rc**2.
    Es= Ec+(Einf-Ec)*10.**logE
    Er= 2.*(Es-Ec) 
    vR= numpy.sqrt(x*Er)
    vz= numpy.sqrt((1-x)*Er)
    o= Orbit([Rc,vR,(10**logL)/Rc,0.,vz,0.])
    return o

def estimate_delta(orbit, ts= numpy.linspace(0.,20.,1001), pot=MWPotential2014):
    '''
    estimate delta for a given orbit instance
    '''
    orbit.integrate(ts,pot,method='symplec4_c')
    return estimateDeltaStaeckel(pot, orbit.R(ts),orbit.z(ts))

def azimuthal_period(orbit, pot=MWPotential2014, delta=0.375):
    ''' 
    find the azimuthal period of an orbit using a given estimate of delta
    '''
    aAS= actionAngleStaeckel(pot=pot, delta=delta)
    _,_,_,_,Op,_= aAS.actionsFreqs(orbit) # computes actions/freqs, only need azimuthal
    if Op == 9999.99:
        return np.nan
    Tp= 2.*numpy.pi/Op
    return Tp
    
def estimate_orbit_params(orbit, pot=MWPotential2014, delta=0.375, c=True):
    '''
    use the staeckel approximation to estimate orbital parameters for an orbit
    '''
    try:
        if not c:
            aASS = actionAngleStaeckelSingle(orbit,pot=pot,delta=delta)
            umin, umax= aASS.calcUminUmax()
            vmin= aASS.calcVmin()
            rperi= bovy_coords.uv_to_Rz(umin,numpy.pi/2.,delta=aASS._delta)[0]
            rap_tmp, zmax = bovy_coords.uv_to_Rz(umax,vmin,delta=aASS._delta)
            rap = numpy.sqrt(rap_tmp**2.+zmax**2.)
            e = (rap-rperi)/(rap+rperi)
            return [rperi, rap, zmax, e]
        else:
            aAS= actionAngleStaeckel(pot=pot,delta=delta)
            e, zmax, rperi, rap = aAS.EccZmaxRperiRap(orbit)
            return [rperi, rap, zmax, e]
    except UnboundError:   
        return [np.nan, np.nan, np.nan, np.nan]
    
def params_along_orbit(orbit, ts, pot=MWPotential2014, delta=0.375):
    '''
    estimate the orbit parameters along an orbit
    '''
    orbit.integrate(ts,pot,method='symplec4_c')
    Rs = orbit.R(ts)
    vRs = orbit.vR(ts)
    vTs = orbit.vT(ts)
    zs = orbit.z(ts)
    vzs = orbit.vz(ts)
    phis = orbit.phi(ts)
    aAS= actionAngleStaeckel(pot=pot,delta=delta) # for example
    es, zms, rps, ras = aAS.EccZmaxRperiRap(Rs,vRs,vTs,zs,vzs,phis)
    return es, zms, rps, ras

