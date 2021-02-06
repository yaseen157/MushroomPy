# -*- coding: utf-8 -*-

"""
MushroomPy Module

Generate plots of detonation velocities from data as produced by the following program(s):
VTF > AMROC > Clawpack > Euler Chemistry > 1D > Detonation Tube

Updated January 2021
Tested with:
    Python 3.8, Windows 10

Author Email: bdmc1n17@soton.ac.uk
"""

__author__ = "Bruno del Mazo Canaleta"


import m1d_amroc as m1d
import numpy as np
from matplotlib import pyplot as plt
from sdtoolbox.postshock import CJspeed, PostShock_fr, PostShock_eq
import os
import glob
import warnings


# Only creates the raw plots of the prescribed cases and parameters.
def rawplots(study, caselist, mechanismlist, hdrlist):
    for case in caselist:
        for mech in mechanismlist:
            for hdr in hdrlist:
                study.export_rawdataplots(case, mech, hdr)

# Produces 
def framerateSpeed(study, caselist, mechanismlist):
    # declare dictionary to hold the data
    dataDict = {}
    # declare lists with the x over time data for each case
    for case in caselist:
        dataDict[case] = {}
        for mech in mechanismlist:
            finalXposition = list(study.volatilecatalogue[case][mech]["dat_0.txt"]["DataframeObj"]["x"])[-1]
            timelist = study.data_statistics()[case][mech]["ReactionFront"]["time (s)"]
            xlist = study.data_statistics()[case][mech]["ReactionFront"]["x"]
            # calculate a list of velocities until the reaction front reaches the end
            # also shorten the list of x values for plotting purposes to fit the velocity data
            velocitylist = []
            newxlist = []
            delta = 1 # number of "frames" used to calculate velocity
            i = delta
            while i < len(xlist)-1:
                if xlist[i]<finalXposition or xlist[i+1]<finalXposition:
                    velocity = (xlist[i]-xlist[i-delta]) / (timelist[i]-timelist[i-delta]) /100 # divide by 100 because x in cm
                    velocitylist.append(velocity)
                    newxlist.append(xlist[i])
                i += 1
            dataDict[case][mech] = np.array([newxlist, velocitylist])
    return dataDict

# Plots the "framerate" wave velocities
def frameratePlotter(study, caselist, mechanismlist):
    dataDict = framerateSpeed(study, caselist, mechanismlist)
    for case in caselist:
        for mech in mechanismlist:
            fig, ax = plt.subplots()
            ax.plot(dataDict[case][mech][0],dataDict[case][mech][1])
            ax.set_xlabel("Distance [cm]")
            ax.set_ylabel("Detonation wave velocity [m/s]")
            ax.grid()
            fig.savefig(f'C:/Users/serio/Desktop/AMROC DATA/YaseensCode/{case}_{mech}.png', bbox_inches='tight')

# density method but taking max near-post-discontinuity value for rho and u
def densityMethodMax(study, caselist, mechanismlist):

    velocityDict = {}
    for case in caselist:
        velocityDict[case] = {}
        for mech in mechanismlist:
            # generate list of datastep.key strings (e.g. dat_0.txt) and the list of positions of the reaction front
            xlist = study.data_statistics()[case][mech]["ReactionFront"]["x"]
            finalXposition = list(study.volatilecatalogue[case][mech]["dat_0.txt"]["DataframeObj"]["x"])[-1]
            dataSteps = [stepstring for stepstring in study.volatilecatalogue[case][mech].keys()]
            
            # shorten the list of useful data points
            newxlist = []
            newsteplist = []
            i = 3   # which data step to start from, excluding as many as needed for comparison across discontinuity
            while i < len(xlist)-1:
                if xlist[i]<finalXposition or xlist[i+1]<finalXposition:
                    newxlist.append(xlist[i])
                    newsteplist.append(dataSteps[i])
                i += 1
            
            # generate lists of densities for continuity equation ahead and behind discontinuity
            # as well as gas velocities behind the discontinuity
            ahead_densitylist = []
            behind_densitylist = []
            gasvelocity = []
            waveVelocity = []
            for i in range(len(newsteplist)):
                allXpositions = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["x"]
                allDensities = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Density"]
                allVelocityUs = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Velocityu"]
                for j in range(1,len(allXpositions)):
                    if allXpositions[j] == newxlist[i]:
                        ahead_densitylist.append(allDensities[j])
                        # find max value of rho and u from near behind shock
                        temp_dens = []; temp_u = []
                        k = 1
                        while allXpositions[j-k] >= max(newxlist[i]-3, 0.5):    # newxlist[i]-3 means it searches for a max up to 3 cm behind the shock, improves accuracy
                            temp_dens.append(allDensities[j-k])
                            temp_u.append(allVelocityUs[j-k])
                            k += 1
                        behind_densitylist.append(max(temp_dens))
                        gasvelocity.append(max(temp_u))

                # calculate the wave speed at each point from continuity (density) eq
                waveV = behind_densitylist[i] * gasvelocity[i] / (behind_densitylist[i] - ahead_densitylist[i])
                waveVelocity.append(waveV)
            # put together as array in a dictionary
            velocityDict[case][mech] = np.array([newxlist, waveVelocity])
    
    return velocityDict

# calculating wave velocity using the eq for the Rayleigh line
def rayleighMethodMax(study, caselist, mechanismlist):

    velocityDict = {}
    for case in caselist:
        velocityDict[case] = {}
        for mech in mechanismlist:
            # generate list of datastep.key strings (e.g. dat_0.txt) and the list of positions of the reaction front
            xlist = study.data_statistics()[case][mech]["ReactionFront"]["x"]
            finalXposition = list(study.volatilecatalogue[case][mech]["dat_0.txt"]["DataframeObj"]["x"])[-1]
            dataSteps = [stepstring for stepstring in study.volatilecatalogue[case][mech].keys()]
            
            # shorten the list of useful data points
            newxlist = []
            newsteplist = []
            i = 3   # which data step to start from, excluding as many as needed for comparison across discontinuity
            while i < len(xlist)-1:
                if xlist[i]<finalXposition or xlist[i+1]<finalXposition:
                    newxlist.append(xlist[i])
                    newsteplist.append(dataSteps[i])
                i += 1
            
            # generate lists of densities and pressures for continuity equation ahead and behind discontinuity
            # as well as gas velocities behind the discontinuity
            ahead_densitylist = []; ahead_pressurelist = []
            behind_densitylist = []; behind_pressurelist = []
            waveVelocity = []
            for i in range(len(newsteplist)):
                allXpositions = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["x"]
                allDensities = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Density"]
                allPressures = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Pressure"]
                for j in range(len(allXpositions)):
                    if allXpositions[j] == newxlist[i]:
                        ahead_densitylist.append(allDensities[j])
                        ahead_pressurelist.append(allPressures[j])
                        # find max value of rho and p from 3 cm behind shock
                        temp_dens = []; temp_pres = []
                        k = 1
                        while allXpositions[j-k] > max(newxlist[i]-3, 0.5):
                            temp_dens.append(allDensities[j-k])
                            temp_pres.append(allPressures[j-k])
                            k += 1
                        behind_densitylist.append(max(temp_dens))
                        behind_pressurelist.append(max(temp_pres))
                # calculate wave velocity with relevant parameters
                waveV = ( (behind_pressurelist[i] - ahead_pressurelist[i]) / (1/ahead_densitylist[i] - 1/behind_densitylist[i]) / ahead_densitylist[i]**2 )**0.5
                waveVelocity.append(waveV)
            # put together as an array in the dictionary
            velocityDict[case][mech] = np.array([newxlist, waveVelocity])
    return velocityDict

# Calculating the theoretical CJ speed using the SDToolbox
def CJmethod(study, caselist, mechanismlist):
    velocityDict = {}
    for case in caselist:
        velocityDict[case] = {}
        for mechanism in mechanismlist:
            finalXposition = list(study.volatilecatalogue[case][mechanism]["dat_0.txt"]["DataframeObj"]["x"])[-1]
            _, p2, moles, _ = case.split('_')   # pressures in kpa
            dil, _, _ = mechanism.split('_')
            T1 = 298    # kelvin
            q = 'C2H4:1. O2:3 ' + dil + ':' + moles
            mech = 'gri30_highT.cti'
            cj = CJspeed(P1=int(p2)*1000, T1=T1, q=q, mech=mech, fullOutput=False)
            xvalues = np.linspace(0,finalXposition,10)
            cjvalues = np.zeros(10)
            for i in range(len(cjvalues)):
                cjvalues[i] = cj
            velocityDict[case][mechanism] = np.array([xvalues, cjvalues])
    return velocityDict


# comparison between different methods in graphical form
def compare(study, caselist, mechanismlist):
    # Check that the output directory exists, create it if not
    detVelocity_path = os.path.join(study.output_path, "DetonationVelocity")
    if not os.path.exists(detVelocity_path):
        warnmsg = f"Path does not exist, generating a new 'DetonationVelocity' directory: ({detVelocity_path})"
        warnings.warn(warnmsg)
        os.makedirs(detVelocity_path)
    for mechanism in mechanismlist:
        mech_path = os.path.join(detVelocity_path, mechanism)
        if not os.path.exists(mech_path):
            warnmsg = f"Path does not exist, generating a new mechanism directory: ({mech_path})"
            warnings.warn(warnmsg)
            os.makedirs(mech_path)
    # theoretical CJ speed taken as reference
    cjDict = CJmethod(study,caselist, mechanismlist)
    # framerate method
    frameDict = framerateSpeed(study, caselist, mechanismlist)
    # density method selecting the max near-post-disc values
    densityDictMax = densityMethodMax(study, caselist, mechanismlist)
    # Rayleigh line method selecting the max near-post-disc values
    rayleightDictMax = rayleighMethodMax(study, caselist, mechanismlist)
    # plotting it all together
    for case in caselist:
        for mech in mechanismlist:
            # Calculate Rayleigh method error at x=300cm
            target = 300
            lst = list(rayleightDictMax[case][mech][0])
            closest = lst[min(range(len(lst)), key = lambda i: abs(lst[i]-target))]
            k = lst.index(closest)
            rayleighV = rayleightDictMax[case][mech][1][k]
            error_percent = round( abs(cjDict[case][mech][1][0] - rayleighV) / cjDict[case][mech][1][0] * 100, 2)
            # plotting
            fig, ax = plt.subplots()
            ax.plot(cjDict[case][mech][0],cjDict[case][mech][1], label="CJ speed", color='red')
            ax.plot(frameDict[case][mech][0],frameDict[case][mech][1], label="framerate", color='black')
            ax.plot(densityDictMax[case][mech][0],densityDictMax[case][mech][1], linestyle='-.', label="density", color='green')
            ax.plot(rayleightDictMax[case][mech][0],rayleightDictMax[case][mech][1], linestyle='-.', label="Rayleigh", color='blue')
            ax.set_xlabel("Distance [cm]")
            ax.set_ylabel("Detonation wave velocity [m/s]")
            ax.set_title(f"Rayleigh: {error_percent}% error at x=300 cm")
            framerate_half = frameDict[case][mech][1][round(len(frameDict[case][mech][1])/2):] # list of values for the second half of the tube
            ymax = max( max(framerate_half), max(cjDict[case][mech][1]) )*1.1
            ymin = min(framerate_half)*0.75
            ax.set_ylim([ymin, ymax])
            ax.grid()
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            fig.savefig(f'{detVelocity_path}/{mech}/{case}_{mech}.png', bbox_inches='tight')
            print(f"{case} plotted for {mech}.")


# This function finds all the available raw data in a prescribed directory and
# runs the compare() function on all of it
def AllIn(): 
    study = m1d.DetonationTube()

    # get all the data folder names from the 1D_AMROC_data folder
    pathcaselist = glob.glob(f'{study.simdata_path}/*')
    caselist = [os.path.basename(i) for i in pathcaselist]
    #initialise mushroomPy
    study = m1d.DetonationTube()

    for case in caselist:
        # get all the data for the relevant case
        study.case_read(pressurestudy = case)
        # get a list of the available mechanism folders inside the particular case folder
        pathmechanismlist = glob.glob(f'{study.simdata_path}/{case}/*')
        mechanismlist = [os.path.basename(i) for i in pathmechanismlist]
        # run compare() to plot all mechanisms for a given case
        compare(study, [case], mechanismlist) # second argument must be a list, even though here we only pass a single case




"""
To use this code you need to define a list of "cases", i.e. 'X_Y_Z_data' where
X and Y are driver and driven pressures (in kpa) and Z is moles of diluent.

Must also define a list of "mechanisms", e.g. "Ar_GRI_red2", which will specify 
which diluent it is and which chemistry code's output to use.

Then read through every case in the list and use whatever function you need.

for example:

study1 = m1d.DetonationTube()
caselist = ["50_10_4_data"]
mechanismlist = ["Ar_GRI_red2"]
for case in caselist:
    study1.case_read(pressurestudy = case)

compare(study1, caselist, mechanismlist)
"""



#study1 = m1d.DetonationTube()
""" caselist = ["10_10_2_data", "10_10_4_data", "10_10_6_data", "10_10_8_data", "10_10_10_data",
            "20_20_2_data", "20_20_4_data", "20_20_8_data", "20_20_10_data"]  """
#caselist = ["100_100_1_data"]

#mechanismlist = ["Ar_GRI_red2"]
#mechanismlist = ["Ar_GRI_red2", "Ar_C2H4_Jachi", "N2_GRI_red2", "N2_C2H4_Jachi"]

""" for case in caselist:
    study1.case_read(pressurestudy = case) """


#compare(study1, caselist, mechanismlist)



""" if __name__ == "__main__":
    AllIn()
 """