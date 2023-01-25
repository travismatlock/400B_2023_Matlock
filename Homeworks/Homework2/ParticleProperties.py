# Author: Travis Matlock
# ASTR 400B 1/24/23
# This script obtains the magnitudes of the distance and velocity for a certain
# particle as well as its mass. It also attached units to these quantities.


# Import modules and functions
import numpy as np
import astropy.units as u
import math
from ReadFile import Read

# This function takes a file name, particle type, and particle number (starting
# from 1). The file name is a string, particularly 'MW_nnn.txt'. The particle
# type is not used yet. The particle number is an integer and starts with the
# first particle. Index corrections are built into the code. It returns the
# magnitude of the distance from the origin in kpc, the magnitude of the
# velocity in km/s, and the mass in solar masses of the specified particle.
def ParticleInfo(filename, type, num):
    # Grab variables from ReadFile. total and time are not yet used.
    time, total, data = Read(filename)
    # Calculate the magnitude of the distance vector
    distance = math.hypot(data['x'][num-1], data['y'][num-1], data['z'][num-1])
    # Assign units and round to three decimal places
    distance = np.around(distance * u.kpc, 3)
    # Calculate the magnitude of the velocity vector
    speed = math.hypot(data['vx'][num-1], data['vy'][num-1], data['vz'][num-1])
    # Assign units and round to three decimal paces
    speed = np.around(speed * u.km / u.s, 3)
    # Extract the mass of and assign units to the particle
    mass = data['m'][num-1] * u.solMass
    # Return finished distance, speed, and mass
    return distance, speed, mass


distance, speed, mass = ParticleInfo('MW_000.txt', 'NA', 100)
print("3D distance: "+str(distance))
print("3D velocity: "+str(speed))
print("Mass: "+str(mass))
print("3D distance in light years: "+str(distance.to(u.lyr)))