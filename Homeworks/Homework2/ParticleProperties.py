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
# type is an integer 1, 2, or 3. The particle number is an integer and starts 
# with the first particle. Index corrections are built in the code. It returns 
# the magnitude of the distance from the origin in kpc, the magnitude of the
# velocity in km/s, and the mass in solar masses of the specified particle.
def ParticleInfo(filename, ptype, num):
    # Grab variables from ReadFile. total and time are not yet used.
    time, total, data = Read(filename)
    # Create a filter that selects particles only with the specified type.
    temp = np.where(data['type'] == ptype)
    # Create an array based on data with only particles of the same type
    type_sorted = data[temp]    
    # Calculate the magnitude of the distance vector
    distance = math.hypot(type_sorted['x'][num-1], type_sorted['y'][num-1], type_sorted['z'][num-1])
    # Assign units and round to three decimal places
    distance = np.around(distance * u.kpc, 3)
    # Calculate the magnitude of the velocity vector
    speed = math.hypot(type_sorted['vx'][num-1], type_sorted['vy'][num-1], type_sorted['vz'][num-1])
    # Assign units and round to three decimal paces
    speed = np.around(speed * u.km / u.s, 3)
    # Extract the mass of and assign units to the particle
    mass = type_sorted['m'][num-1] * u.solMass
    # Return finished distance, speed, and mass
    return distance, speed, mass