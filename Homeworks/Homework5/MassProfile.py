# Author: Travis Matlock
# ASTR 400B 2/16/2023

# Import Relevant modules, functions, constants
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G
# Convert G to necessary units
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

class MassProfile:
# Class to define mass profiles, rotation curves, and Hernquist profiles
    def __init__(self, galaxy, snap):
        '''
        This function initializes the class in order to obtain mass profiles
        and rotation curves of the galaxies

        Parameters
        ----------
        galaxy : string
            Name of the galaxy
        snap : int
            The SnapNumber (time) of the simulation

        Returns
        -------
        None.

        '''
        # Convert snap input to a 3-digit snapnumber string preceded by 0s
        ilbl = '000'+str(snap)
        ilbl = ilbl[-3:]
        # Convert inputs to .txt file name format
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        # Calculate time, number of particles, and data arrays of .txt file
        self.time, self.total, self.data = Read(self.filename)
        # Define methods for calculating x, y, z positions and mass
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']
        # Store galaxy name
        self.gname = galaxy
        
    def MassEnclosed(self, ptype, radii):
        '''
        This function calculates the mass profile of a component of a galaxy
        
        Parameters
        ----------
        ptype : int
            Component to calculat (1 = halo, 2 = disk, 3 = bulge)
        radii : array
            Each radius to calculate enclosed mass for
            
        Returns
        -------
        mass_profile : array
            An Array of mass enclosed w/ elements corresponding to radii elements
            
        '''
        # Call CenterOfMass and initialize
        center_of_mass = CenterOfMass(self.filename, ptype)
        # Calculate COM position within .1 kpc
        com_position = center_of_mass.COM_P(0.1)
        # Create empty array to store masses
        mass_profile = np.array([])
        # Create a method to calculate the distance from a particle to COM
        self.r = np.sqrt((self.x-com_position[0])**2 +
                              (self.y-com_position[1])**2 +
                              (self.z-com_position[2])**2)
        # Loop through each radius in radii
        for r in radii:
            # Create a filter for all particles within radius
            self.r_index = np.where(self.r < (r*u.kpc))
            # Sort particles by this filter
            self.r_sorted = self.data[self.r_index]
            # Create a filter for all particles of this component
            self.t_index = np.where(self.r_sorted['type'] == ptype)
            # Store an array of masses filtered by ptype within looped r
            masses = self.r_sorted['m'][self.t_index]
            # Calculate the total mass of all these masses
            total_mass = np.sum(masses)
            # Add this total mass to the mass profile array
            mass_profile = np.append(mass_profile, total_mass)
        # Return mass profile array, converting to Msun.
        return mass_profile * 1e10 * u.Msun
    
    def MassEnclosedTotal(self, radii):
        '''
        This function calculates the mass profile of all components of a galaxy
        
        Parameters
        ----------
            radii: array
                Each radius to calculate enclosed mass for
            
        Returns
        -------
            total_mass: Array
                The total mass profile array of the galaxy
                
        '''
        # Calculate halo and disk elements
        halo = self.MassEnclosed(1,radii)
        disk = self.MassEnclosed(2,radii)
        # M33 has no bulge component, so it needs an exception
        if self.gname == 'M33':
            # Set bulge to an array of 0s
            bulge = np.zeros(len(radii))
        else:
            # Calculate bulge element
            bulge = self.MassEnclosed(3,radii)
        # Add all component mass profiles
        total_mass = np.add(halo, disk)
        total_mass = np.add(total_mass, bulge)
        # Return mass profile of all components
        return total_mass
    
    def HernquistMass(self, r, a, Mhalo):
        '''
        This function calculates the Hernquist mass of a galaxy
        
        Parameters
        ----------
            r: float
                The radius to calculate M at
            a: float
                The scale radius of the galaxy
            Mhalo: float
                The mass of halo components of the galaxy
                
        Returns
        -------
            M: float
                Hernquist mass of the galaxy
                
        '''
        # Calculate and return the Hernquist mass
        M = Mhalo*r**2 / (a+r)**2
        return M
    
    def CircularVelocity(self, ptype, radii):
        '''
        This function calculates the circular velocity of a component at a
        specified radius
        
        Parameters
        ----------
            ptype: int
                Specific galactic component
            radii: array
                An array of different radii to calculate circular velocity at
                
        Returns
        -------
            v_c: float
                The circular velocities of the galaxy at a specified distance
        
        '''
        # Calculate the mass profile of this specific component
        mass_arr = self.MassEnclosed(ptype, radii)
        # Calculate circular velocity of this component at each radius
        v_c = np.sqrt(G * mass_arr / (radii*u.kpc))
        return v_c
    
    def CircularVelocityTotal(self, radii):
        '''
        This function calculates the circular velocity of all mass enclosed
        at given radii

        Parameters
        ----------
        radii : array
            The array of radius values to calculate v_c for

        Returns
        -------
        v_c : array
            The array of circular velocities at each radius

        '''
        # Calculate a total mass array for each radius
        mass_arr = self.MassEnclosedTotal(radii)
        # Calculate the circular velocities of all components at each radius
        v_c_total = np.sqrt(G*mass_arr/(r*u.kpc))
        return v_c_total

    
    def HernquistVCirc(self, radii, a, Mhalo):
        '''
        This function calculates the circular velocity at different radii
        based on their enclosed Hernquist mass
        
        Parameters
        ----------
            radii: array
                The radii to calculate M at
            a: float
                The scale radius of the galaxy
            Mhalo: float
                The mass of halo components of the galaxy
                
        Returns 
        -------
            v_c_hernquist: array
                An array of circular velocities at each radius
                
        '''
        # Calculate the Hernquist mass array at different radii
        mass_arr = self.HernquistMass(radii, a, Mhalo)
        # Calculate the Hernquist circular velocity
        v_c_hernquist = np.sqrt(G * mass_arr / (radii*u.kpc))
        # Round to 2 decimal places and return
        return np.round(v_c_hernquist, 2)
 
# Create an array of radii from .1 kpc to 30 kpc in .1 kpc intervals
r = np.arange(.1,30.1,.1)
# Initialize class for this galaxy name and snap
MW = MassProfile('MW', 0)
# Calculate the halo, disk, bulge, Hernquist, and total masses
MW_halo = MW.MassEnclosed(1,r)
MW_disk = MW.MassEnclosed(2,r)
MW_bulge = MW.MassEnclosed(3,r)
MW_total = MW.MassEnclosedTotal(r)
MW_hernquist = MW.HernquistMass(r,28.3,1.975E12)
# Plot the mass profiles for components and totals
fig, ax = plt.subplots()
# Use log scale for y axis
plt.semilogy()
# Create labels
ax.set(title='Milky Way Mass Profile', xlabel='Radius [kpc]', ylabel='Encompassed Mass [Msun]')
ax.plot(r, MW_halo, label='Halo')
ax.plot(r, MW_disk, 'purple', label='Disk')
ax.plot(r, MW_bulge, 'g', label='Bulge')
ax.plot(r, MW_total, 'black', label='Total')
ax.plot(r, MW_hernquist, 'r:', label='Hernquist (a=28.3)' )
ax.legend()

# Initialize class for this galaxy name and snap
M31 = MassProfile('M31', 0)
# Calculate the halo, disk, bulge, Hernquist, and total masses
M31_halo = M31.MassEnclosed(1,r)
M31_disk = M31.MassEnclosed(2,r)
M31_bulge = M31.MassEnclosed(3,r)
M31_total = M31.MassEnclosedTotal(r)
M31_hernquist = M31.HernquistMass(r,28.8,1.9721E12)
# Plot the mass profiles for components and totals
fig, ax = plt.subplots()
# Use log scale for y axis
plt.semilogy()
# Create labels
ax.set(title='M31 Mass Profile', xlabel='Radius [kpc]', ylabel='Encompassed Mass [Msun]')
ax.plot(r, M31_halo, label='Halo')
ax.plot(r, M31_disk, 'purple', label='Disk')
ax.plot(r, M31_bulge, 'g', label='Bulge')
ax.plot(r, M31_total, 'black', label='Total')
ax.plot(r, M31_hernquist, 'r:', label='Hernquist (a=28.8)' )
ax.legend()

# Initialize class for this galaxy name and snap
M33 = MassProfile('M33', 0)
# Calculate the halo, disk, bulge, Hernquist, and total masses
M33_halo = M33.MassEnclosed(1,r)
M33_disk = M33.MassEnclosed(2,r)
M33_total = M33.MassEnclosedTotal(r)
M33_hernquist = M33.HernquistMass(r,15.8,1.87E11)
# Plot the mass profiles for components and totals
fig, ax = plt.subplots()
# Use log scale for the y axis
plt.semilogy()
# Create labels
ax.set(title='M33 Mass Profile', xlabel='Radius [kpc]', ylabel='Encompassed Mass [Msun]')
ax.plot(r, M33_halo, label='Halo')
ax.plot(r, M33_disk, 'purple', label='Disk')
ax.plot(r, M33_total, 'black', label='Total')
ax.plot(r, M33_hernquist, 'r:', label='Hernquist (a=15.8)' )
ax.legend()



# Calculate the component and total circular velocities
MW_halo = MW.CircularVelocity(1,r)
MW_disk = MW.CircularVelocity(2,r)
MW_bulge = MW.CircularVelocity(3,r)
MW_total = MW.CircularVelocityTotal(r)
MW_hernquist = MW.HernquistVCirc(r,28.3,1.975E12)
# Plot the rotation curves for components and totals
fig, ax = plt.subplots()
# Use log scale for y axis
plt.semilogy()
# Create labels
ax.set(title='Milky Way Rotation Curve', xlabel='Radius [kpc]', ylabel='Circular velocity [km/s]')
ax.plot(r, MW_halo, label='Halo')
ax.plot(r, MW_disk, 'purple', label='Disk')
ax.plot(r, MW_bulge, 'g', label='Bulge')
ax.plot(r, MW_total, 'black', label='Total')
ax.plot(r, MW_hernquist, 'r:', label='Hernquist (a=28.3)' )
ax.legend()

# Initialize class for this galaxy name and snap
M31 = MassProfile('M31', 0)
# Calculate the halo, disk, bulge, Hernquist, and total circular velocities
M31_halo = M31.CircularVelocity(1,r)
M31_disk = M31.CircularVelocity(2,r)
M31_bulge = M31.CircularVelocity(3,r)
M31_total = M31.CircularVelocityTotal(r)
M31_hernquist = M31.HernquistVCirc(r,28.8,1.921E12)
# Plot the rotation curves for components and totals
fig, ax = plt.subplots()
# Use log scale for y axis
plt.semilogy()
# Create labels
ax.set(title='M31 Rotation Curve', xlabel='Radius [kpc]', ylabel='Circular Velocity [km/s]')
ax.plot(r, M31_halo, label='Halo')
ax.plot(r, M31_disk, 'purple', label='Disk')
ax.plot(r, M31_bulge, 'g', label='Bulge')
ax.plot(r, M31_total, 'black', label='Total')
ax.plot(r, M31_hernquist, 'r:', label='Hernquist (a=28.8)' )
ax.legend()

# Initialize class for this galaxy name and snap
M33 = MassProfile('M33', 0)
# Calculate the halo, disk, bulge, Hernquist, and total circular velocities
M33_halo = M33.MassEnclosed(1,r)
M33_disk = M33.MassEnclosed(2,r)
M33_total = M33.MassEnclosedTotal(r)
M33_hernquist = M33.HernquistVCirc(r,15.8,1.87E11)
# Plot the rotation curves for components and totals
fig, ax = plt.subplots()
# Use log scale for the y axis
plt.semilogy()
# Create labels
ax.set(title='M33 Rotation Curve', xlabel='Radius [kpc]', ylabel='Circular Velocity [km/s]')
ax.plot(r, M33_halo, label='Halo')
ax.plot(r, M33_disk, 'purple', label='Disk')
ax.plot(r, M33_total, 'black', label='Total')
ax.plot(r, M33_hernquist, 'r:', label='Hernquist (a=15.8)')
ax.legend()