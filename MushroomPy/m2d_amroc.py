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

import os
import warnings

import cantera as ct
from sdtoolbox.postshock import CJspeed, PostShock_fr
from sdtoolbox.znd import zndsolve 

import mpy_config as mpyconf


class prep_2d:
    """
    **TWO-dimensional** simulation preparation tool. Intended to be used to automatically generate solver.in and init.dat files
    for prescribed conditions for use with the 2D computations (StatDet, MovingDet or Obstacle).
    """
    def __init__(self, P1=10000., T1=298., diluent='Ar', dilpercent=50.):

        # Find the dataroot directory, or prompt user to choose one if non-existent
        self.data_path = mpyconf.parse_configfile()["dataroot"]
        # Python module name
        self.module = os.path.basename(__file__).replace(".py", "")
        # Python class name
        self.classname = type(self).__name__

        # Define input and output directories
        self.output_path = os.path.join(self.data_path, "_output", self.module, self.classname)
        self.mushroompath = os.path.dirname(os.path.realpath(__file__))
        self.input_path = os.path.join(self.mushroompath, "templates")

        # Check that templates folder exists
        if not os.path.exists(self.input_path):
            warnmsg = f"Templates path does not exist, generating a new one: {self.input_path}\n\nPlease, provide init.dat and solver.in example files in the specified directory.\n"
            warnings.warn(warnmsg)
            os.makedirs(self.input_path)
        
        # Check that the templates folder actually has templates
        """
        WRITE LATER
        """

        # Check that prep output folder exists
        if not os.path.exists(self.output_path):
            warnmsg = f"Path for {self.classname} folder does not exist, generating a new one: {self.output_path}"
        
        # Use initial conditions to find CJ speed 
        self.P1 = P1
        self.T1 = T1
        self.diluent = diluent
        self.dilpercent = dilpercent
        self.dilmoles = self.perc2mol(self.dilpercent)
        self.SDTmech = 'gri30_highT.cti'
        self.q = 'C2H4:1. O2:3 ' + self.diluent + ':' + str(self.dilmoles)
        self.cj_speed = CJspeed(self.P1, self.T1, self.q, self.SDTmech, fullOutput=False)

    # Convert percentage dilution into molar fractions
    def perc2mol(self, percent, fuelmols=1, oxmols=3):
        dilmoles = percent * (fuelmols+oxmols) / (100 - percent)
        return dilmoles

    # Calculate induction length for a given case
    def induction_length(self, cj_speed=False, P1=False, T1=False, q=False, mech=False, t_end=1e-5):
        # Assign default values innit
        if cj_speed is False:
            cj_speed = self.cj_speed
        if P1 is False:
            P1 = self.P1
        if T1 is False:
            T1 = self.T1
        if q is False:
            q = self.q
        if mech is False:
            mech = self.SDTmech
        # Set up gas object
        gas1 = ct.Solution(mech)
        gas1.TPX = T1,P1,q
        # Find post-shock conditions
        gas = PostShock_fr(cj_speed, P1, T1, q, mech)
        # Solve ZND ODEs to find the width of the ZND "plateau"
        znd_out = zndsolve(gas,gas1,cj_speed,t_end=t_end,advanced_output=True)
        plateau_length = znd_out['ind_len_ZND']
        
        return plateau_length



    """ # Returns the minimum number of cells required in the base mesh to
    # achieve at least k cells between the shock front and the combustion front.
    def minCells(self, P1=10000, T1=298, diluent='Ar', dilpercent=50, length=40, k=5, maxLevel=3, RF=2):
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
        return cells """

    # Function that generates an init file with desired conditions
    def generateInitFile(self, labFrame=True):

        # get lines from the template file
        initf = open(f"{self.input_path}/init.dat", "r")
        lines = initf.readlines()
        # modify the elements of the first line which define initial conditions
        initConditions = lines[0].split(" ")
        initConditions[0] = f"{round(self.T1,1)}d0"
        initConditions[1] = f"{round(self.P1*1e-3,2)}d0"
        separator = ' '
        lines[0] = separator.join(initConditions)
        # modify elements of second line, defining frame and wave speed
        cjline = lines[1].split(' ')
        cjline[1] = f"{round(self.cj_speed,2)}d0" # wave velocity 
        if labFrame:
            cjline[2] = '0'
        else:
            cjline[2] = '1'
        lines[1] = separator.join(cjline)
        # modify elements of third line, defining initial moral fractions of different species
        molefracs = lines[2].split(' ')
        if self.diluent == 'Ar':
            molefracs[10] = f"{round(self.dilmoles,2)}d0"
        elif self.diluent == 'N2':
            molefracs[9] = f"{round(self.dilmoles,2)}d0"
        else:
            warnmsg = "Unknown diluent"
            warnings.warn(warnmsg)
        lines[2] = separator.join(molefracs)

        # write new init file in a folder for the corresponding conditions
        separator = '_'
        """ case_directory = separator.join([[self.diluent, ]]) """

        print(f"{self.output_path}\\init.dat")

        



"""data_path = mpyconf.parse_configfile()["dataroot"]"""





"""MIGHT HAVE TO ADD A BODGE FACTOR OF *3 TO DELTA X"""
#minCells(dilpercent=0)