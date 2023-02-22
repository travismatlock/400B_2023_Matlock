# Author: Travis Matlock
# ASTR 400B 2/21/23

# Import modules and function
import numpy as np
import astropy
import astropy.units as u
import matplotlib.pyplot as plt
import sys
from ReadFile import Read
from CenterOfMass import CenterOfMass

def OrbitCOM(galaxy, start=0, end=800, n=5):
    '''
    This function generates .txt files that describe the position and velocity
    of the specified galaxy's COM for a specified duration of the simulation'
    
    INPUTS
    ------
        galaxy: 'str'
            The desired galaxy for which to track COM
        start: 'int'
            The snap number at which to start tracking COM. Default is 0
        end: 'int'
            The snap number at which to stop tracking COM. Default is 800,
            which is the end of the simulation.
        n: 'int'
            Step size by which to increase snap number each iteration
            
    OUTPUTS
    -------
        fileout: 'Orbit_[GALAXY].txt'
            A .txt file detailing time, component velocities, and component
            positions of the COM of the galaxy, saved to wdir.
    '''
    # Name the final .txt file
    fileout = 'Orbit_'+galaxy+'.txt'
    # Set acceptable error tolerance for COM calculation
    delta = .1
    # Set factor by which COM volume will be reduced during calculation
    volDec = 2
    # M33 needs a higher factor due to tidal stripping
    if galaxy == 'M33':
        volDec = 4
    # Create an array of every snap number that will contribute to orbit tracking
    snap_ids = np.arange(start, end+1, n)
    # Store amount of snap numbers that are being tracked
    snap_dur = np.size(snap_ids)
    # If only 1 snap or fewer is being tracked, exit the code
    if snap_dur <= 1:
        sys.exit()
    # Create an array with 7 columns and as many rows as tracked snap numbers
    # This array will store time, position, and velocity of the COM
    orbit = np.zeros([snap_dur, 7])
    # Loop over each desired snap number. Use indexes for loop count and snap number
    for i, snap_id in enumerate(snap_ids):
        # Uncomment if you desire to track how many times the loops has been run
        # print(i)
        # Format snap number as a 3-digit number
        ilbl = '000'+str(snap_id)
        ilbl = ilbl[-3:]
        # Create the file name of the galaxy and snap number to calculate COM for
        filename = "%s_"%(galaxy) + ilbl + '.txt'
        # Create a CenterOfMass class object for the current snap number.
        # Use disk stars 
        COM = CenterOfMass(filename, ptype=2)
        # Calculate position of COM
        mc_p = (COM.COM_P(delta, volDec))
        # Calculate velocity using COM position components
        mc_v = (COM.COM_V(mc_p[0], mc_p[1], mc_p[2]))
        # Store time of each snap in 1st column. Units are Gyr, but object
        # is dimensionless
        orbit[i,0] = (COM.time/1000).value
        # Store x, y, z of COM in columns 2, 3, 4. Units are kpc but object
        # is dimensionless
        orbit[i,1] = mc_p[0].value
        orbit[i,2] = mc_p[1].value
        orbit[i,3] = mc_p[2].value
        # Store vx, vy, vz of COM in columns 5, 6, 7. Units are km/s but object
        # is dimensionless
        orbit[i,4] = mc_v[0].value
        orbit[i,5] = mc_v[1].value
        orbit[i,6] = mc_v[2].value
    # Format and save the text file to the working directory
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:10s}{:11s}{:11s}{:11s}{:11s}{:11s}{:11s}"
               .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

# Uncomment to generate Orbit_[GALAXY].txt files for specified galaxies
# OrbitCOM('MW')
# OrbitCOM('M31')
# OrbitCOM('M33')

def RelativeProperties(gname1, gname2):
    '''
    This function takes two galaxy names and calculates the relative distance
    and velocity of the two galaxies' CsOM.
    
    INPUTS
    ------
        gname1: 'str'
            First galaxy in set. This one is slightly favored by the code, as
            it is used to calculate time and length, although these values
            are expected be the same.
        gname2: 'str'
            Galaxy to copmare to first galaxy
            
    RETURNS
    -------
        time_data: 'array'
            A 1D array of dimnesionless time values in Gyr.
        sep_data: 'array'
            A 1D array of dimensionless relative distances in kpc.
        rel_vel_data: 'array'
            A 1D array of dimensionless relative velocities in km/s
    '''
    # Open COM tracking files generated from OrbitCOM in read mode
    file1 = open('Orbit_'+gname1+'.txt', 'r')
    file2 = open('Orbit_'+gname2+'.txt', 'r')
    # From files, generate a 2D array where each row represents a snap number,
    # Each column holds specified orbital data.
    file1_data = np.genfromtxt(file1,skip_header=1)
    file2_data = np.genfromtxt(file2,skip_header=1)
    # Create empty arrays to store incremented values
    time_data = np.array([])
    sep_data = np.array([])
    rel_vel_data = np.array([])
    # Loop over each row in Orbit_[GALAXY].txt
    for i in range(len(file1_data)):
        # Calculate the component relative distance values of galaxies' COMs
        x_diff = abs(file1_data[i][1] - file2_data[i][1])
        y_diff = abs(file1_data[i][2] - file2_data[i][2])
        z_diff = abs(file1_data[i][3] - file2_data[i][3])
        # Calculate magnitude of relative distance between COMs
        sep = np.sqrt(x_diff**2+y_diff**2+z_diff**2)
        # Calculate the component relative velocity values of galaxies' COMs
        vx_diff = abs(file1_data[i][4] - file2_data[i][4])
        vy_diff = abs(file1_data[i][5] - file2_data[i][5])
        vz_diff = abs(file1_data[i][6] - file2_data[i][6])
        # Calculate magnitude of relative velocity between COMs
        rel_vel = np.sqrt(vx_diff**2+vy_diff**2+vz_diff**2)
        # Store time, separation, and relative velocity for current snap
        time_data = np.append(time_data, file1_data[i][0])
        sep_data = np.append(sep_data, sep)
        rel_vel_data = np.append(rel_vel_data, rel_vel)
    # Return finished arrays
    return time_data, sep_data, rel_vel_data

# Calculate desired data arrays for the M31-M33 system and for the MW-M31 system
sat_t, sat_sep, sat_rel_vel = RelativeProperties('M31', 'M33')
titans_t, titans_sep, titans_rel_vel = RelativeProperties('MW', 'M31')
# Create a figure with 2 plots and format it according to desires
fig, ax = plt.subplots(1, 2)
fig.tight_layout(h_pad=5)
fig.set_figwidth(10)
fig.set_figheight(5)
#ax[0].semilogy()
# Create overlaid plots of separation and relative velocity for both systems
ax[0].set(title='Relative Distance of Galaxies', xlabel='Time [Gyr]', ylabel='Separation [kpc]')
ax[1].set(title='Relative Velocities of Galaxies', xlabel='Time [Gyr]', ylabel='Relative Velocity [km/s]')
ax[0].plot(sat_t, sat_sep, 'r', label='M31-M33')
ax[0].plot(titans_t, titans_sep, label='MW-M31')
ax[1].plot(sat_t, sat_rel_vel, 'r', label='M31-M33')
ax[1].plot(titans_t, titans_rel_vel, label='MW-M31')
ax[0].legend()
ax[1].legend()