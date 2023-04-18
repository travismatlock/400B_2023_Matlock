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
        # Create COM-frame pos, vel, component and vector arrays for all ptype particles
        self.x = self.x_raw - self.com_pos_vec[0]
        self.y = self.y_raw - self.com_pos_vec[1]
        self.z = self.z_raw - self.com_pos_vec[2]
        self.vx = self.vx_raw - self.com_vel_vec[0]
        self.vy = self.vy_raw - self.com_vel_vec[1]
        self.vz = self.vz_raw - self.com_vel_vec[2]
        self.pos_vec = np.array([self.x,self.y,self.z]).T
        self.vel_vec = np.array([self.vx,self.vy,self.vz]).T
        
    
        
    def RotateFrame(self,posI,velI):
        """a function that will rotate the position and velocity vectors
        so that the disk angular momentum is aligned with z axis. 
        
        PARAMETERS
        ----------
            posI : `array of floats`
                 3D array of positions (x,y,z)
            velI : `array of floats`
                 3D array of velocities (vx,vy,vz)
                 
        RETURNS
        -------
            pos: `array of floats`
                rotated 3D array of positions (x,y,z) such that disk is in the XY plane
            vel: `array of floats`
                rotated 3D array of velocities (vx,vy,vz) such that disk angular momentum vector
                is in the +z direction 
        """
        
        # compute the angular momentum
        L = np.sum(np.cross(posI,velI), axis=0)
        # normalize the vector
        L_norm = L/np.sqrt(np.sum(L**2))
    
    
        # Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
        
        # z unit vector
        z_norm = np.array([0, 0, 1])
        
        # cross product between L and z
        vv = np.cross(L_norm, z_norm)
        s = np.sqrt(np.sum(vv**2))
        
        # dot product between L and z 
        c = np.dot(L_norm, z_norm)
        
        # rotation matrix
        I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
        R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2
    
        # Rotate coordinate system
        pos = np.dot(R, posI.T).T
        vel = np.dot(R, velI.T).T
        
        return pos, vel    
            
    def SunlikeIndices(self):
        '''
        This function defines an annular disk of stars with sunlike orbital radii
        and velocities. It will be run for the initial snapnumber.
        
        OUTPUTS
        -------
        indices: `array of ints`
            Indices of snap file that represent stellar candidates
        '''
        # Create empty storage array
        indices = np.array([])
        # Establish type-sorted array
        temp = self.data[self.type_index]
        # Rotate vectors so that angular momentum is in z direction
        pos_rotated, vel_rotated = self.RotateFrame(self.pos_vec, self.vel_vec)
        # Obtain the radial distance in rotated frame
        xy_dist = np.sqrt(pos_rotated[:,0]**2 + pos_rotated[:,1]**2)
        # Create filter for radial distance
        sunlike_annulus = np.where(np.isclose(xy_dist, 8.178, atol=1.2267))
        # Filter for height from disk
        sunlike_heights = np.where(np.abs(pos_rotated[:,2]) < 0.5)
        # Calculate velocity magnitude and filter
        speed = np.sqrt(vel_rotated[:,0]**2+vel_rotated[:,1]**2+vel_rotated[:,2]**2)
        sunlike_speed = np.where(np.isclose(speed, 242.6, atol=18.96))
        # Create storage array
        indices = np.array([])
        # Loop over all 3 filters
        for i in sunlike_annulus[0]:
            # If an index (of temp) is in all three filters, record it
            if i in sunlike_heights[0] and i in sunlike_speed[0]:
                    indices = np.append(indices, i)
        # Return array of sunlike indices
        return indices
    
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
        speeds = np.sqrt(self.vel_vec[indices,0]**2 + self.vel_vec[indices,1]**2 + self.vel_vec[indices,2]**2)
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
        distances = np.sqrt(self.pos_vec[indices,0]**2 + self.pos_vec[indices,1]**2 + self.pos_vec[indices,2]**2)
        return distances
    
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
        mags = self.SunlikeDistances(indices)
        mags = np.array([mags, mags, mags]).T
        r_hats = self.pos_vec[indices]/mags
        v_rads = np.array([])
        velocities = self.vel_vec[indices]
        for j in range(r_hats.shape[0]):
            v_rad = np.dot(r_hats[j],velocities[j])
            v_rads = np.append(v_rads, v_rad)
        return v_rads
    # Do not need GetProperties as function will be methods for Snap objects.
        
# Create snap object for 'M31_000.txt' and define stellar candidates
# Check these indices at 'M31_800.txt' and create histograms of separation,
# magnitude and radial velocities.
# Loop over snapnumbers in multiples of 5. Find average position and velocity
# vectors at these times and graph as a function of time.
snap_i = Snap('M31_000.txt')
stellar_candidates = snap_i.SunlikeIndices().astype(int)
#print(stellar_candidates)
(snap_i.RadialVelocity(stellar_candidates))

snap_f = Snap('M31_800.txt')
final_orbital_distances = snap_f.SunlikeDistances(stellar_candidates)
fig, ax = plt.subplots()
ax.set(title='Final Orbital distances of stellar candidate stars',
       xlabel='Orbital Distance [kpc]', ylabel='Number of Stellar Candidates')
ax.hist(final_orbital_distances, bins=100, range=(0,60))
plt.vlines(8.178,0,50, 'r')

count = 0
for dist in final_orbital_distances:
    if dist >= 8.178:
        count += 1
print(count/final_orbital_distances.size)

# Obtain time evolution of velocities and position.
snap_ids = np.arange(0,800,20)
# Initialize storage arrays
radial_velocities = np.array([])
speeds = np.array([])
distances = np.array([])
times = np.array([])
for i, snap_id in enumerate(snap_ids):
    # Format snap number as a 3-digit number
    ilbl = '000'+str(snap_id)
    ilbl = ilbl[-3:]
    # Create the file name of the galaxy and snap number to calculate COM for
    filename = 'M31_' + str(ilbl) + '.txt'
    # Create snap object
    current_snap = Snap(filename)
    mean_rad_vel = np.mean(current_snap.RadialVelocity(stellar_candidates))
    mean_speed = np.mean(current_snap.SunlikeSpeeds(stellar_candidates))
    mean_distance = np.mean(current_snap.SunlikeDistances(stellar_candidates))
    radial_velocities = np.append(radial_velocities, mean_rad_vel)
    speeds = np.append(speeds, mean_speed)
    distances = np.append(distances, mean_distance)
    times = np.append(times, current_snap.time.value)

fig, ax = plt.subplots()
ax.plot(times, radial_velocities, 'g')
ax.plot(times, speeds, 'b')
ax.plot(times, distances, 'r')
    
    
    