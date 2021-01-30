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

# Density method
def densityMethod(postdiscposition, study, caselist):
    
    
    # REMOVE FIRST ARGUMENT AFTER THE BEST WAY IS DECIDED

    velocityDict = {}
    for case in caselist:
        # generate list of datastep.key strings (e.g. dat_0.txt) and the list of positions of the reaction front
        xlist = study.data_statistics()[case]["Ar_GRI_red2"]["ReactionFront"]["x"]
        dataSteps = [stepstring for stepstring in study.volatilecatalogue[case]["Ar_GRI_red2"].keys()]
        
        # shorten the list of useful data points
        newxlist = []
        newsteplist = []
        i = 2   # which data step to start from, excluding as many as needed for comparison across discontinuity
        while xlist[i]<349.5 or xlist[i+1]<349.5:
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
            allXpositions = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["x"]
            allDensities = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["Density"]
            allVelocityUs = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["Velocityu"]
            for j in range(len(allXpositions)):
                if allXpositions[j] == newxlist[i]:
                    ahead_densitylist.append(allDensities[j])
                    behind_densitylist.append(allDensities[j-postdiscposition])    # might have to make it i-2 or 3 to account for potential vN
                    gasvelocity.append(allVelocityUs[j-postdiscposition])
            # calculate the wave speed at each point from continuity (density) eq
            waveV = behind_densitylist[i] * gasvelocity[i] / (behind_densitylist[i] - ahead_densitylist[i])
            waveVelocity.append(waveV)
        # put together as array in a dictionary
        velocityDict[case] = np.array([newxlist, waveVelocity])
    
    return velocityDict

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

# Momentum method
def momentumMethod(postdiscposition, study, caselist):
    
    
    # REMOVE FIRST ARGUMENT AFTER THE BEST WAY IS DECIDED

    velocityDict = {}
    for case in caselist:
        # generate list of datastep.key strings (e.g. dat_0.txt) and the list of positions of the reaction front
        xlist = study.data_statistics()[case]["Ar_GRI_red2"]["ReactionFront"]["x"]
        dataSteps = [stepstring for stepstring in study.volatilecatalogue[case]["Ar_GRI_red2"].keys()]
        
        # shorten the list of useful data points
        newxlist = []
        newsteplist = []
        i = 2   # which data step to start from, excluding as many as needed for comparison across discontinuity
        while xlist[i]<349.5 or xlist[i+1]<349.5:
            newxlist.append(xlist[i])
            newsteplist.append(dataSteps[i])
            i += 1
        
        # generate lists of densitiesand pressures for continuity equation ahead and behind discontinuity
        # as well as gas velocities behind the discontinuity
        ahead_densitylist = []; ahead_pressurelist = []
        behind_densitylist = []; behind_pressurelist = []
        gasvelocity = []
        waveVelocity = []
        for i in range(len(newsteplist)):
            allXpositions = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["x"]
            allDensities = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["Density"]
            allVelocityUs = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["Velocityu"]
            allPressures = study.volatilecatalogue[case]["Ar_GRI_red2"][newsteplist[i]]["DataframeObj"]["Pressure"]
            for j in range(len(allXpositions)):
                if allXpositions[j] == newxlist[i]:
                    ahead_densitylist.append(allDensities[j])
                    ahead_pressurelist.append(allPressures[j])
                    behind_densitylist.append(allDensities[j-postdiscposition])    # might have to make it i-2 or 3 to account for potential vN
                    behind_pressurelist.append(allPressures[j-postdiscposition])
                    gasvelocity.append(allVelocityUs[j-postdiscposition])
            # calculate the wave speed at each point from momentum eq
            # eq turns out to be a quadratic, so it's simpler to group terms
            a = ahead_densitylist[i] - behind_densitylist[i]
            b = 2 * behind_densitylist[i] * gasvelocity[i]
            c = ahead_pressurelist[i] - behind_pressurelist[i] - behind_densitylist[i] * gasvelocity[i]**2
            sol1 = (-1*b + (b**2 - 4*a*c)**0.5) / (2*a)
            sol2 = (-1*b - (b**2 - 4*a*c)**0.5) / (2*a)
            print (a, b, c, b**2 -4*a*c, sol1, sol2)

# Momentum method but taking max near-post-discontinuity value for rho and u
# NEEDS WORK, CURRENTLY IT BASICALLY IS THE SAME AS DENSITY METHOD
def momentumMethodMax(study, caselist, mechanismlist):

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
            gasvelocity = []
            sol1list = []; sol2list = []
            for i in range(len(newsteplist)):
                allXpositions = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["x"]
                allDensities = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Density"]
                allVelocityUs = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Velocityu"]
                allPressures = study.volatilecatalogue[case][mech][newsteplist[i]]["DataframeObj"]["Pressure"]
                for j in range(1, len(allXpositions)):
                    if allXpositions[j] == newxlist[i]:
                        ahead_densitylist.append(allDensities[j])
                        ahead_pressurelist.append(allPressures[j])
                        # find max value of rho, p and u from 3 cm behind shock
                        temp_dens = []; temp_u = []; temp_pres = []
                        k = 1
                        while allXpositions[j-k] > max(newxlist[i]-3, 0.5):
                            temp_dens.append(allDensities[j-k])
                            temp_u.append(allVelocityUs[j-k])
                            temp_pres.append(allPressures[j-k])
                            k += 1
                        behind_densitylist.append(max(temp_dens))
                        behind_pressurelist.append(max(temp_pres))
                        gasvelocity.append(max(temp_u))
                # calculate the wave speed at each point from momentum eq
                # eq turns out to be a quadratic, so it's simpler to group terms
                a = ahead_densitylist[i] - behind_densitylist[i]
                b = 2 * behind_densitylist[i] * gasvelocity[i]
                c = ahead_pressurelist[i] - behind_pressurelist[i] - behind_densitylist[i] * gasvelocity[i]**2
                if (b**2 - 4*a*c) < 0:
                    sol1 =  (-1*b) / (2*a)
                    sol2 =  (-1*b) / (2*a)
                else:
                    sol1 = (-1*b + (b**2 - 4*a*c)**0.5) / (2*a)
                    sol2 = (-1*b - (b**2 - 4*a*c)**0.5) / (2*a)
                sol1list.append(sol1)
                sol2list.append(sol2)
            velocityDict[case][mech] = np.array([newxlist, sol1list, sol2list])
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
            p1, p2, moles, leftover = case.split('_')   # pressures in kpa
            dil, leftover1, leftover2 = mechanism.split('_')
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
    # density method with varying delta-x across shock
    #densityDict1 = densityMethod(1, study, caselist, mechanismlist)
    #densityDict2 = densityMethod(2, study, caselist, mechanismlist)
    #densityDict3 = densityMethod(3, study, caselist, mechanismlist)
    # density method selecting the max near-post-disc values
    densityDictMax = densityMethodMax(study, caselist, mechanismlist)
    # Momentum method selecting the max near-post-disc values
    momentumDictMax = momentumMethodMax(study, caselist, mechanismlist)
    # Rayleigh line method selecting the max near-post-disc values
    rayleightDictMax = rayleighMethodMax(study, caselist, mechanismlist)
    # plotting it all together
    for case in caselist:
        for mech in mechanismlist:
            fig, ax = plt.subplots()
            ax.plot(cjDict[case][mech][0],cjDict[case][mech][1], label="CJ speed", color='red')
            ax.plot(frameDict[case][mech][0],frameDict[case][mech][1], label="framerate", color='black')
            #ax.plot(densityDict1[case][mech][0],densityDict1[case][mech][1], linestyle='--', label="density (i-1)")
            #ax.plot(densityDict2[case][mech][0],densityDict2[case][mech][1], linestyle='--', label="density (i-2)")
            #ax.plot(densityDict3[case][mech][0],densityDict3[case][mech][1], linestyle='--', label="density (i-3)")
            ax.plot(densityDictMax[case][mech][0],densityDictMax[case][mech][1], linestyle='-.', label="density MAX", color='green')
            #ax.plot(momentumDictMax[case][mech][0],momentumDictMax[case][mech][1], linestyle='-.', label="momentum MAX sol1", color='orange')
            #ax.plot(momentumDictMax[case][mech][0],momentumDictMax[case][mech][2], linestyle='-.', label="momentum MAX sol2", color='purple')
            ax.plot(rayleightDictMax[case][mech][0],rayleightDictMax[case][mech][1], linestyle='-.', label="Rayleigh MAX", color='blue')
            ax.set_xlabel("Distance [cm]")
            ax.set_ylabel("Detonation wave velocity [m/s]")
            ax.set_ylim([0, max( max(frameDict[case][mech][1]), max(cjDict[case][mech][1]) )*1.1])
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

compare(study1, caselist)
"""



#study1 = m1d.DetonationTube()
#caselist = ["100_100_0_data", "100_100_2_data", "100_100_4_data"]
#caselist = ["50_10_2_data", "50_10_4_data", "50_10_6_data", "50_10_18_data"]
""" caselist = ["10_10_2_data", "10_10_4_data", "10_10_6_data", "10_10_8_data", "10_10_10_data",
            "20_20_2_data", "20_20_4_data", "20_20_8_data", "20_20_10_data"]  """

#caselist = ["20_20_4_data", "20_20_6_data"]

#mechanismlist = ["Ar_GRI_red2"]
#mechanismlist = ["N2_C2H4_Jachi"]
#mechanismlist = ["Ar_GRI_red2", "Ar_C2H4_Jachi", "N2_GRI_red2", "N2_C2H4_Jachi"]

""" for case in caselist:
    study1.case_read(pressurestudy = case)
 """

#hdrlist = ["Velocityu", "Pressure", "SpeedofSound"]
#hdrlist = ["Density"]
#rawplots(study1, caselist, mechanismlist, hdrlist)


#frameratePlotter(study1, caselist, mechanismlist)

#compare(study1, caselist, mechanismlist)

#momentumMethodMax(study1, caselist, mechanismlist)
#rayleighMethodMax(study1, caselist, mechanismlist)
#CJmethod(study1, caselist, mechanismlist)

""" for case in caselist:
    for mech in mechanismlist:
        xlist = study1.data_statistics()[case][mech]["ReactionFront"]["x"]
        dataSteps = [stepstring for stepstring in study1.volatilecatalogue[case][mech].keys()]
        newsteplist = []
        newxlist = []
        i = 2   # which data step to start from, excluding as many as needed for comparison across discontinuity
        while xlist[i]<349.5 or xlist[i+1]<349.5:
            newxlist.append(xlist[i])
            newsteplist.append(dataSteps[i])
            i += 1
        print(len(newxlist))
        allXpositions = study1.volatilecatalogue[case][mech][newsteplist[23]]["DataframeObj"]["x"]
        for j in range(len(allXpositions)):
            if allXpositions[j] == newxlist[23]:
                #print(allXpositions[j-3], allXpositions[j-2], allXpositions[j-1], allXpositions[j])
                print(allXpositions[j] - allXpositions[j-3]) 
 """

AllIn()
