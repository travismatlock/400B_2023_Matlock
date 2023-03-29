# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 

# Completed by: Travis Matlock
# ASTR 400B 3/19/23


# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# import the RelativeProperties to obtain alternatively simulated data.
from OrbitCOM import RelativeProperties


# # M33AnalyticOrbit

class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        #        """ **** ADD COMMENTS """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.fileout = filename[0:3]+'_LF.txt'
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass('M33_000.txt', 2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        posM33 = (M33_COM).COM_P(.1, 4)
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        velM33 = (M33_COM).COM_V(posM33[0],posM33[1],posM33[2])
        
        posM33 = posM33.value
        velM33 = velM33.value
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31_COM = CenterOfMass('M31_000.txt', 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        posM31 = (M31_COM).COM_P(.1, 2)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        velM31 = (M31_COM).COM_V(posM31[0],posM31[1],posM31[2])
        posM31 = posM31.value
        velM31 = velM31.value
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r = posM33 - posM31
        self.v = velM33 - velM31
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 #kpc
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass('M31_000.txt', 2).value  * 1e12# Msun
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1 #kpc
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass('M31_000.txt', 3).value *1e12 # Msun
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 #kpc
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass('M31_000.txt', 1).value *1e12 #Msun
    
    
    def HernquistAccel(self, M, r_a, r): 
        '''
        This function calculates the Hernquist profile acceleration on a galaxy.
        
        INPUTS
        ------
            M: `float`
                Total mass of component of galaxy in Msun
            r_a: `float`
                Hernquist profile scale radius in kpc
            r: `array`
                x,y,z position of COM of orbiting galaxy in kpc
                
        RETURNS
        -------
            Hern: `array`
                Acceleration vector based on a Hernquist profile
        '''
        
        ### **** Store the magnitude of the position vector
        rmag = np.linalg.norm(r)
        
        ### *** Store the Acceleration
        
        Coef = -self.G*M / (rmag * (r_a + rmag)**2)
        Hern = Coef*r 
        return Hern
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        '''
        This function calculate the Miyamoto-Nagai profile acceleration on a galaxy
        
        INPUTS
        ------
            M: `float`
                Mass of component of galaxy in Msun
            r_d: `float`
                Scale radius of disk in kpc
            r: `array`
                x,y,z positions of COM of galaxy in kpc
                
        RETURNS
        -------
            a: `array`
                Acceleration vector based on the Miyamoto-Nagai disk profile.
        '''
        
        ### Acceleration **** follow the formula in the HW instructions
        # Calculate parameters
        z_d = self.rdisk/5.0
        R = np.sqrt(r[0]**2+r[1]**2)
        B = r_d + np.sqrt(r[2]**2+z_d**2)
        # Calculate acceleration components
        a_x = -self.G*M*r[0] / (R**2+B**2)**1.5
        a_y = -self.G*M*r[1] / (R**2+B**2)**1.5
        a_z = -self.G*M*B*r[2] / (R**2+B**2)**1.5 / (r[2]**2+z_d**2)
        # Return acceleration vector
        a = np.array([a_x,a_y,a_z])
       
        return a
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r):
        '''
        This function calculates the total acceleration on M31.
        
        INPUTS
        ------
            r: `array`
                x,y,z positions of the COM of the galaxy
            
        RETURNS
        -------
            a: `array`
                Acceleration vector on the COM based on all particle types
                
        '''
        # Calculate acceleration vector for each galactic component
        a_halo = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        a_bulge = self.HernquistAccel(self.Mbulge, self.rbulge, r)
        a_disk = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
        # Calculate acceleration vector for each vector component
        a_x = a_halo[0] + a_bulge[0] + a_disk[0]
        a_y = a_halo[1] + a_bulge[1] + a_disk[1]
        a_z = a_halo[2] + a_bulge[2] + a_disk[2]
        # Return the total acceleration vector
        a = np.array([a_x, a_y, a_z])
        return a
    
    
    
    def LeapFrog(self, dt, r, v):
        '''
        This function integrates the equation of motion over a differential time.
        
        INPUTS
        ------
            dt: `float`
                Timestep in Gyr
            r: `array`
                Current x,y,z positions of COM
            v: `array`
                Current velocity vector of COM
        
        RETURNS
        -------
            rnew: `array`
                Integrated position vector of COM in dt
            vnew: `array`
                Integrated velocity vector of COM in dt
        '''
        # predict the position at the next half timestep
        rhalf = np.array([r[0]+v[0]*dt/2, r[1]+v[1]*dt/2, r[2]+v[2]*dt/2])
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        a = self.M31Accel(rhalf)
        vnew = np.array([v[0]+a[0]*dt, v[1]+a[1]*dt, v[2]+a[2]*dt])
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = np.array([rhalf[0]+vnew[0]*dt/2, rhalf[1]+vnew[1]*dt/2, rhalf[2]+vnew[2]*dt/2])
        
        return rnew, vnew
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        '''
        This function integrates the orbit of the galaxy over a time interval in differential timesteps

        Parameters
        ----------
        t0 : `float`
            Initial time (Gyr)
        dt : `float`
            Timestep (Gyr)
        tmax : `float`
            Final time (Gyr)

        Returns
        -------
        None.

        '''

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r), *tuple(self.v)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t < tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            r = orbit[i-1,1:4]
            #print(r)
            v = orbit[i-1,4:8]
            #print(v)
            rnew, vnew = self.LeapFrog(dt, r, v)
            #print(rnew, vnew)
            #print('\n')
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i, 1:4] = rnew
            orbit[i, 4:7] = vnew
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            i += 1
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
        
        
        # write the data to a file
        np.savetxt(self.fileout, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
# Initialize object of class M33AnalyticOrbit
M33 = M33AnalyticOrbit('M33_000.txt')
# Integrate from 0 to 10 Gyr
M33.OrbitIntegration(0, 0.001, 10)
# Obtain integrated data
M33Orbit = np.genfromtxt('M33_LF.txt', names=True)
# Initialize plot
fig, ax = plt.subplots(1,2)
fig.tight_layout(w_pad=5)
# Calculate total separation magnitude
integrated_sep = np.sqrt(M33Orbit['x']**2+M33Orbit['y']**2+M33Orbit['z']**2)
# Obtain desired values from alternative simulation
sat_t, sat_sep, sat_rel_vel = RelativeProperties('M31', 'M33')
# Graph separations
ax[0].set(title='Separation of M31-M33', xlabel='time [Gyr]', ylabel='Separation [kpc]')
ax[0].plot(sat_t, sat_sep, 'r', label='Defined Orbit')
ax[0].plot(M33Orbit['t'][:-1], integrated_sep[:-1], label='Analytical Orbit')
ax[0].legend()
# Graph relative velocities
ax[1].set(title='Relative Velocities of M31-M33', xlabel='Time [Gyr]', ylabel='Relative Velocity [kpc/Gyr]')
ax[1].plot(sat_t, sat_rel_vel, 'r', label='Defined Velocity')
rel_vels = np.sqrt(M33Orbit['vx']**2+M33Orbit['vy']**2+M33Orbit['vz']**2)
ax[1].plot(M33Orbit['t'][:-1], rel_vels[:-1], label='Analytical Velocity')
ax[1].legend(loc='upper left')