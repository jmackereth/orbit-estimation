import numpy as np
from galpy.potential import SCFPotential, DiskSCFPotential

class mySCFPotential(SCFPotential):
    def _R2deriv(self,R,z,phi=0.,t=0.):
        dR= 1e-8
        return (self._Rforce(R,z) - self._Rforce(R+dR,z))/dR
    
    def _z2deriv(self,R,z,phi=0.,t=0.):
        dz = 1e-8
        return (self._zforce(R,z) - self._zforce(R,z+dz))/dz
    
    def _Rzderiv(self,R,z,phi=0.,t=0.):
        dR = 1e-8
        return (self._zforce(R,z) - self._zforce(R+dR,z))/dR
        
        
class myDiskSCFPotential(DiskSCFPotential):
    def _R2deriv(self,R,z,phi=0.,t=0.):
        dR= 1e-8
        return (self._Rforce(R,z) - self._Rforce(R+dR,z))/dR
    
    def _z2deriv(self,R,z,phi=0.,t=0.):
        dz = 1e-8
        return (self._zforce(R,z) - self._zforce(R,z+dz))/dz
    
    def _Rzderiv(self,R,z,phi=0.,t=0.):
        dR = 1e-8
        return (self._zforce(R,z) - self._zforce(R+dR,z))/dR
        