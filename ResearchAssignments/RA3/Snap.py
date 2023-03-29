# Author: Travis Matlock
# ASTR 400B 3/29/23

''' This is a draft using a class for snapshot data. Doing this will allow
me to easily access all particle data at different snapshots.'''

''' My research topic is the kinematic evolution of sunlike disk stars in M31'''

''' Question: How do the kinematics of Sun analogs change? '''


# Import modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

# Import functions
from ReadFile_Soln import Read
from CenterOfMass_Soln import CenterOfMass

class Snap:
    def __init__(self, filename, ptype=2):
        '''
        This function initializes an object of class Snap

        Parameters
        ----------
        filename : `str`
            Snapshot file to initialize as Snap object
        ptype : `int`, optional
            Integer representing which particle type to examine. 
            The default is 2 for disk particles
        '''
        print('init '+filename)
        # Read in file data
        self.time, self.total, self.data = Read(filename)
        # Create an index to sort by specified ptype
        self.type_index = np.where(self.data['type'] == ptype)
        # Create methods for particle data similar to CenterOfMass
        self.x_raw = self.data['x'][self.type_index]
        self.y_raw = self.data['y'][self.type_index]
        self.z_raw = self.data['z'][self.type_index]
        self.vx_raw = self.data['vx'][self.type_index]
        self.vy_raw = self.data['vy'][self.type_index]
        self.vz_raw = self.data['vz'][self.type_index]
        # Create COM object and get COM pos, vel vectors
        com = CenterOfMass(filename, ptype)
        com_pos_vec = com.COM_P(0.1)
        self.com_vel_vec = com.COM_V(com_pos_vec[0],com_pos_vec[1], com_pos_vec[2]).value
        self.com_pos_vec = com_pos_vec.value
         # Create pos, vel vector array, type sorted and adjusted to COM frame
        pos_vec = np.array([])
        vel_vec = np.array([])
        # Loop over each type selected particle
        for i in range(len(self.x_raw)):
            # Create vectors in MW (appx) frame
            curr_pos_vec = np.array([self.x_raw[i], self.y_raw[i], self.z_raw[i]])
            curr_vel_vec = np.array([self.vx_raw[i], self.vy_raw[i], self.vz_raw[i]])
            # Amend these vectors into COM frame
            curr_pos_vec -= self.com_pos_vec
            curr_vel_vec -= self.com_vel_vec
            # Append COM frame vectors to array of all COM pos, vel vectors.
            pos_vec = np.append(pos_vec, curr_pos_vec)
            vel_vec = np.append(vel_vec, curr_vel_vec)
        # Create methods for COM frame pos, vel vector array for all particles
        self.pos_vec = pos_vec.reshape((len(self.x_raw),3))
        self.vel_vec = vel_vec.reshape((len(self.vx_raw),3))
    
    def SunlikeIndices(self):
        print('indices')
        '''
        This function defines an annulus of stars with sunlike orbital radii
        and velocities. It will be run for the initial snapnumber.
        
        OUTPUTS
        -------
        indices: `array of ints`
            Indices of snap file that represent stellar candidates
        '''
        # Create an empty storage array
        indices = np.array([])
        # Loop over particles
        j = 0
        for i in range(self.pos_vec.shape[0]):
            # Calculate magnitude of pos, vel vectors
            dist = np.linalg.norm(self.pos_vec[i])
            speed = np.linalg.norm(self.vel_vec[i])
            # If these match stellar properties, store their index
            if 7.5 < dist and dist < 8.5 and 220 < speed and speed < 250:
                indices = np.append(indices, i)
        # Return list of sunlike indices.
        return indices
    
    def RadialVelocity(self, indices):
        '''
        This function computes the radial velocity array for indexed particles only.
        
        INPUTS
        ------
        indices: `array of ints`
            Indices of snap file that represent stellar candidates
            
        OUTPUTS
        -------
        v_rads: `array of floats`
            Radial velocity components for each particle
        '''
        # Create empty storage array
        v_rads = np.array([])
        # Loop over selected stars
        for i in indices:
            # Calculate magnitude of current r vector
            mag = np.linalg.norm(self.pos_vec[i])
            # Calculate r hat
            r_hat = self.pos_vec[i] / mag
            # Dot r hat, velcotiy vector to get radial velocity
            v_rad = np.dot(r_hat,self.vel_vec[i])
            # Store radial velocity
            v_rads = np.append(v_rads, v_rad)
        # Return 1d array of each indexed star's v_rad.
        return v_rads
    
    def SunlikeSpeeds(self, indices):
        #print('started speeds')
        '''
        This function calculates the magnitude of the velocity vector for indexed
        particles only
        
        INPUTS
        ------
        indices: `array of ints`
            Indices of snap file that represent stellar candidates
            
        OUTPUTS
        -------
        speeds: `array of floats`
            Magnitude of velocity for each stellar candidate
        '''
        speeds = np.array([])
        for i in indices:
            speed = np.linalg.norm(self.vel_vec[int(i)])
            speeds = np.append(speeds, speed)
        return speeds
    
    def SunlikeDistances(self, indices):
        #print('started distances')
        '''
        This function calculates the magnitude of the position vector for indexed
        particles only
        
        INPUTS
        ------
        indices: `array of ints`
            Indices of snap file that represent stellar candidates
            
        OUTPUTS
        -------
        distances: `array of floats`
            Magnitude of position for each stellar candidate
        '''
        distances = np.array([])
        for i in indices:
            #print(i)
            distance = np.linalg.norm(self.pos_vec[int(i)])
            distances = np.append(distances, distance)
        return distances
        
    # Do not need GetProperties as function will be methods for Snap objects.
        
# Create snap object for 'M31_000.txt' and define stellar candidates
# Check these indices at 'M31_800.txt' and create histograms of separation,
# magnitude and radial velocities.
    'This does of course assume indices for each particle stay the same'
# Loop over snapnumbers in multiples of 5. Find average position and velocity
# vectors at these times and graph as a function of time.
snap_i = Snap('M31_000.txt')
stellar_candidates = snap_i.SunlikeIndices()
snap_f = Snap('M31_800.txt')
final_orbital_distances = snap_f.SunlikeDistances(stellar_candidates)
fig, ax = plt.subplots()
ax.set(title='Final Orbital distances of stellar candidate stars',
       xlabel='Orbital Distance [kpc]', ylabel='Number of Stellar Candidates')
ax.hist(final_orbital_distances, bins=100, range=(0,60))
ax.plot(np.full(25, 8), np.linspace(0,25,25), 'r:', linewidth='1')