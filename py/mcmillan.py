import numpy as np
from galpy.potential import DiskSCFPotential, NFWPotential, \
    SCFPotential, scf_compute_coeffs_axi
from galpy.util import bovy_conversion
from SCF_derivs import mySCFPotential, myDiskSCFPotential

ro = 8.21
vo = 233.1
sigo = bovy_conversion.surfdens_in_msolpc2(vo=vo, ro=ro)
rhoo = bovy_conversion.dens_in_msolpc3(vo=vo, ro=ro)

#gas disk parameters (fixed in McMillan (2017)...)
Rd_HI = 7/ro
Rm_HI = 4/ro
zd_HI = 0.085/ro
Sigma0_HI = 53.1/sigo
Rd_H2 = 1.5/ro
Rm_H2 = 12/ro
zd_H2 = 0.045/ro
Sigma0_H2 = 2180/sigo

#parameters of best-fitting model in McMillan (2017)
#stellar disks
Sigma0_thin = 896/sigo
Rd_thin = 2.5/ro
zd_thin = 0.3/ro
Sigma0_thick = 183/sigo
Rd_thick = 3.02/ro
zd_thick = 0.9/ro
#bulge
rho0_bulge = 98.4/rhoo
r0_bulge = 0.075/ro
rcut = 2.1/ro
#DM halo
rho0_halo = 0.00854/rhoo
rh = 19.6/ro

def gas_dens(R,z, separate=False):
    '''
    gas disk density as a f(R,z). Not sure if this handles R/z == 0 correctly... need to avoid divide by 0.
    '''
    if hasattr(R, "__iter__"):
        HI_dens = np.empty(len(R))
        H2_dens = np.empty(len(R))
        Requal0 = R == 0.
        HI_dens[Requal0] = 0.
        H2_dens[Requal0] = 0.
        HI_dens[~Requal0] = Sigma0_HI/(4*zd_HI)*np.exp(-Rm_HI/R[~Requal0]-R[~Requal0]/Rd_HI)*(sech(z[~Requal0]/(2*zd_HI)))**2
        H2_dens[~Requal0] = Sigma0_H2/(4*zd_H2)*np.exp(-Rm_H2/R[~Requal0]-R[~Requal0]/Rd_H2)*(sech(z[~Requal0]/(2*zd_H2)))**2
        if separate:
            return HI_dens, H2_dens
        return HI_dens+H2_dens   
    else:
        if R == 0.:
            HI_dens = Sigma0_HI/(4*zd_HI)*np.exp(-np.inf-R/Rd_HI)*(sech(z/(2*zd_HI)))**2
            H2_dens = Sigma0_H2/(4*zd_H2)*np.exp(-np.inf-R/Rd_H2)*(sech(z/(2*zd_H2)))**2
            HI_dens = 0.
            H2_dens = 0.
        else:
            HI_dens = Sigma0_HI/(4*zd_HI)*np.exp(-Rm_HI/R-R/Rd_HI)*(sech(z/(2*zd_HI)))**2
            H2_dens = Sigma0_H2/(4*zd_H2)*np.exp(-Rm_H2/R-R/Rd_H2)*(sech(z/(2*zd_H2)))**2
        if separate:
            return HI_dens, H2_dens
        return HI_dens+H2_dens

def stellar_dens(R,z):
    '''
    stellar disk density as a f(R,z)
    '''
    thin_dens = Sigma0_thin/(2*zd_thin)*np.exp(-np.fabs(z)/zd_thin-R/Rd_thin)
    thick_dens = Sigma0_thick/(2*zd_thick)*np.exp(-np.fabs(z)/zd_thick-R/Rd_thick)
    return thin_dens+thick_dens

def bulge_dens(R,z):
    '''
    bulge density as a f(R,z)
    '''
    rdash = np.sqrt(R**2+(z/0.5)**2)
    dens = rho0_bulge/(1+rdash/r0_bulge)**1.8*np.exp(-(rdash/rcut)**2)
    return dens

def NFW_dens(R,z):
    '''
    NFW profile as f(R,z) - not actually used to generate the potential
    '''
    r = np.sqrt(R**2+z**2)
    x = r/rh
    dens = rho0_halo/(x*(1+x)**2)
    return dens

def bulge_gas_dens(R,z):
    '''
    combine bulge, gas and stellar density models
    '''
    return bulge_dens(R,z)+gas_dens(R,z)+stellar_dens(R,z)

def gas_stellar_dens(R,z):
    '''
    combine only gas and stellar density
    '''
    return gas_dens(R,z)+stellar_dens(R,z)

def tot_dens(R,z):
    '''
    total density (including NFW halo) as f(R,z)
    '''
    return gas_dens(R,z)+stellar_dens(R,z)+bulge_dens(R,z)+NFW_dens(R,z)

def sech(x):
    return 1./np.cosh(x)


#dicts used in DiskSCFPotential 
sigmadict = [{'type':'exp','h':Rd_HI,'amp':Sigma0_HI, 'Rhole':Rm_HI},
             {'type':'exp','h':Rd_H2,'amp':Sigma0_H2, 'Rhole':Rm_H2},
             {'type':'exp','h':Rd_thin,'amp':Sigma0_thin, 'Rhole':0.},
             {'type':'exp','h':Rd_thick,'amp':Sigma0_thick, 'Rhole':0.}]

hzdict = [{'type':'sech2', 'h':zd_HI},
          {'type':'sech2', 'h':zd_H2},
          {'type':'exp', 'h':0.3/ro},
          {'type':'exp', 'h':0.9/ro}]

#generate separate disk and halo potential - and combined potential
McMillan_bulge=\
    mySCFPotential(Acos=scf_compute_coeffs_axi(bulge_dens,20,10,a=0.1)[0],
                 a=0.1,ro=ro,vo=vo)
McMillan_disk = myDiskSCFPotential(dens=lambda R,z: gas_stellar_dens(R,z),
                                 Sigma=sigmadict, hz=hzdict,
                                 a=2.5, N=30, L=30,ro=ro,vo=vo)
McMillan_halo = NFWPotential(amp = rho0_halo*(4*np.pi*rh**3),
                             a = rh,ro=ro,vo=vo)
McMillan2017 = [McMillan_disk,McMillan_halo,McMillan_bulge]


                                                            
