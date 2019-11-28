# -*- coding: ascii -*-
"""
Creates a protoplanetary disk around a sun-like star
"""
from __future__ import print_function
import numpy

from matplotlib import pyplot

from amuse.community.fi.interface import Fi

from amuse.units import units
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles


def make_map(sph, N=500, L=1):

    x, y = numpy.indices((N + 1, N + 1))

    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.AU(x)
    y = units.AU(y)
    z = units.AU(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rho.reshape((N + 1, N + 1))

    return numpy.transpose(rho)


if __name__ == "__main__":

    N = 20000
    rentang = numpy.linspace(0, 100, 41)
    for j in rentang:
		
	    tend = j | units.yr
	    umur = str(j)
	    Mstar = 3. | units.MSun

	    convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
	    proto = ProtoPlanetaryDisk(N, convert_nbody=convert, densitypower=1.5, Rmin=4, Rmax=50, q_out=1.)
	    gas = proto.result
	    gas.h_smooth = 0.01 | units.AU

	    sun = Particles(1)
	    sun.mass = Mstar
	    sun.radius = 2. | units.AU
	    sun.x = 0. | units.AU
	    sun.y = 0. | units.AU
	    sun.z = 0. | units.AU
	    sun.vx = 0. | units.kms
	    sun.vy = 0. | units.kms
	    sun.vz = 0. | units.kms

	    sph = Fi(convert)

	    sph.parameters.use_hydro_flag = True
	    sph.parameters.radiation_flag = False
	    sph.parameters.self_gravity_flag = True
	    sph.parameters.gamma = 1.
	    sph.parameters.isothermal_flag = True
	    sph.parameters.integrate_entropy_flag = False
	    sph.parameters.timestep = 0.25 | units.yr
	    
	    sph.gas_particles.add_particles(gas)
	    sph.particles.add_particles(sun)

	    sph.evolve_model(tend)

	    L = 130
	    rho = make_map(sph, N=500, L=L)
	    sph.stop()
	    pyplot.figure(figsize=(8, 8))
	    pyplot.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)), extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15, cmap='inferno')
	    pyplot.title(umur+' tahun')
	    pyplot.colorbar(label='densitas', orientation='vertical', fraction=0.046, pad=0.04,)
	    pyplot.xlabel('au')
	    pyplot.ylabel('au')
	    pyplot.tight_layout()
	    pyplot.savefig('proplyds/proplyd_20000N_3Msun_' +umur+'_tahun.png', dpi = 300)
    #pyplot.show()
