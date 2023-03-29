# Travis Matlock
# ASTR 400B 3/27/23

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

# Import functions
from ReadFile_Soln import Read
from CenterOfMass_Soln import CenterOfMass


def DefineSunLikeParticles(filename, ptype=2):
    '''
    This function determines which indices in a snapshot file represent Sun-like stars.
    
    INPUTS
    ------
    filename: `
    
    OUTPUTS
    -------
    sunlike_indices: `list of int`
        A list of indices at which a particle has sunlike properties at snap 000.
    '''
    # Read in the file
    time, total, data = Read(filename)
    # Create an index to sort by specified ptype
    type_index = np.where(data['type'] == ptype)
    # Create arrays of all position, velocity components
    x = data['x'][type_index]
    y = data['y'][type_index]
    z = data['z'][type_index]
    vx = data['vx'][type_index]
    vy = data['vy'][type_index]
    vz = data['vz'][type_index]
    # Obtain unitless position and velocity vectors of the COM
    COM = CenterOfMass(filename, ptype)
    com_pos_vec = COM.COM_P(0.1)
    com_vel_vec = COM.COM_V(com_pos_vec[0],com_pos_vec[1], com_pos_vec[2]).value
    com_pos_vec = com_pos_vec.value
    # Create empty storage list
    sunlike_indices = []
    # Loop over each particle in type-sorted array
    for i in range(len(x)):
        # Calculate position / velocity vector of every particle.
        pos_vec = np.array([x[i], y[i], z[i]])
        vel_vec = np.array([vx[i], vy[i], vz[i]])
        # Calculate speed and distance of each particle relative to COM
        dist = np.linalg.norm(pos_vec - com_pos_vec)
        speed = np.linalg.norm(vel_vec - com_vel_vec)
        # If kinematic properties fall in range of solar values, record their type sorted index
        # Currently using 8 +/- 0.5 kpc and 235 +/- 15 km/s
        if 7.5 < dist and dist < 8.5 and 220 < speed and speed < 250:
            sunlike_indices.append(i)
        # Return list of all indices with solar properties.
    return sunlike_indices
  

    # Check indexed particles at snap 800
def GetProperties(filename, indices, ptype = 2):
    print(indices)
    time, total, data = Read(filename)
    type_index = np.where(data['type'] == ptype)
    COM = CenterOfMass(filename, ptype)
    com_pos_vec = COM.COM_P(0.1)
    com_vel_vec = COM.COM_V(com_pos_vec[0],com_pos_vec[1], com_pos_vec[2]).value
    com_pos_vec = com_pos_vec.value
    x = data['x'][type_index]
    y = data['y'][type_index]
    z = data['z'][type_index]
    vx = data['vx'][type_index]
    vy = data['vy'][type_index]
    vz = data['vz'][type_index]
    pos_vec = np.array([x,y,z])
    sep = np.array([])
    rel_vel = np.array([])
    for i in indices:
        r = np.array([x[i], y[i], z[i]]) - com_pos_vec
        rv = np.array([vx[i],vy[i],vz[i]]) - com_vel_vec
        sep = np.append(sep, np.linalg.norm(r))
        rel_vel = np.append(rel_vel, np.linalg.norm(rv))
        
    
    
    # Plot frequency histogram of final distance and relative velocity
        # Relative velocity needs to be magnitude and radial
    fig, ax = plt.subplots()
    ax.hist(sep, bins=100, range=(0.,10.))
    # Check indexed particles at every 5th snap number
    
    # Calculate average velocity and distance at each snap number
        # Average is better than median because individual fluctuations are more visible.
    
    # Plot median distance and velocity as functions of time.
        
# Biggest help I need is with organization. Should I make this as a class?
# I tried but decided that it would be too time-inefficient since you'd need a new object
# for every snapnumber. However, if there is a way to create a class that would not need this,
# that would be great. I could create a 3d array of all snap data with time as the third dimension.
# Would this be too big to do anything with or would its capabilites be more worth it than redoing
# these operations for each snap file?

sunlike_indices = DefineSunLikeParticles('M31_000.txt')
GetProperties('M31_800.txt',sunlike_indices)