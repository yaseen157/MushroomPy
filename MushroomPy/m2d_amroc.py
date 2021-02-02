# -*- coding: utf-8 -*-

"""
MushroomPy Module

Updated January 2021
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
Co-Author Email: bdmc1n17@soton.ac.uk
"""

__author__ = "Yaseen Reza"

import cantera as ct
from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve 

import mpy_config as mpyconf

"""Functions used to manipulate the solver and init files"""
# Convert percentage dilution into molar fractions
def perc2mol(percent, fuelmols=1, oxmols=3):
    dilmoles = percent * (fuelmols+oxmols) / (100 - percent)
    return dilmoles

# Returns the minimum number of cells required in the base mesh to
# achieve at least k cells between the shock front and the combustion front.
def minCells(P1=10000, T1=298, diluent='Ar', dilpercent=50, length=40, k=5, maxLevel=3, RF=2):
    # First must simulate the ZND structure to predict the distance between fronts, aka the "plateau"
    mech = 'gri30_highT.cti'
    # Define the mixture in molar quantities
    dilmoles = perc2mol(dilpercent)
    q = 'C2H4:1. O2:3 ' + diluent + ':' + str(dilmoles)
    # Obtain predicted CJ speed
    cj_speed = CJspeed(P1, T1, q, mech, fullOutput=False)
    # Set up gas object
    gas1 = ct.Solution(mech)
    gas1.TPX = T1,P1,q
    # Find post-shock conditions
    gas = PostShock_fr(cj_speed, P1, T1, q, mech)
    # Solve ZND ODEs to find the width of the ZND "plateau"
    znd_out = zndsolve(gas,gas1,cj_speed,t_end=1e-5,advanced_output=True)
    plateau_length = znd_out['ind_len_ZND']
    # Find required deltaX (in cm)
    dX = plateau_length/k *100
    print(dX)
    # Find the number of cells in base mesh to achieve dX given maxLevel
    cells = length/(dX*RF**(maxLevel-1))
    print(cells)
    return cells

# Function that generates an init file with desired conditions
def generateInitFile(P1=10000, T1=298, diluent='Ar', dilpercent=50, labFrame=True):
    # First must calculate the CJ speed for the prescribed conditions
    mech = 'gri30_highT.cti'
    # Define the mixture in molar quantities
    dilmoles = perc2mol(dilpercent)
    q = 'C2H4:1. O2:3 ' + diluent + ':' + str(dilmoles)
    # Obtain predicted CJ speed
    cj_speed = CJspeed(P1, T1, q, mech, fullOutput=False)



"""data_path = mpyconf.parse_configfile()["dataroot"]"""





"""MIGHT HAVE TO ADD A BODGE FACTOR OF *3 TO DELTA X"""
minCells(dilpercent=0)