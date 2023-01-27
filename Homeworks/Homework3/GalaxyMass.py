# Author: Travis Matlock
# ASTR 400B 1/26/2023

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


def ComponentMass(filename, ptype):
    """ This function calculates the total mass of any component of any galaxy.
        
        Inputs:
            filename: 'string'
                The name of the .txt file to be read
            ptype: 'int'
                The type of particle (component) for which to calculate mass
            
        Outputs:
            running_mass: 'float'
                The mass of the specified galactic component, rounded to three
                decimal places and in units of 1e12 solar masses.
    """
    # Grab data from the text file
    time, total, data = Read(filename)
    # Create a filter that selects particles only with the specified type.
    temp = np.where(data['type'] == ptype)
    # Create an array based on data with only particles of the same type
    type_sorted = data[temp]
    # Create a cumulative sum
    running_mass = 0
    # Loop through each particle in the type_sorted array
    for i in range(len(type_sorted)):
        # Add the mass to the cumulative sum of masses
        running_mass += type_sorted['m'][i]
    # Define a unit for 1e12 solar masses
    twelfth_msun = u.def_unit('1e12solMass', 1e12 * u.solMass)
    # Convert the sum to 1e12 solar masses and assign units
    running_mass = running_mass * 10**-2 * twelfth_msun
    # Return the cumulastive sum rounded to three decimal places.
    return np.around(running_mass, 3)