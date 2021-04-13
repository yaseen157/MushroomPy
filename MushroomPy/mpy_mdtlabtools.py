# -*- coding: utf-8 -*-

"""
MushroomPy Module

Tools for laboratory work with the Southampton Modular Detonation Tube.

Updated April 2021
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
"""

__author__ = "Yaseen Reza"


import numpy as np


def calc_networkvolume(pipelengths_dict):
    """Use this function to calculate the volume of a pipe.

    **Parameters:**

    pipelengths_dict
        dictionary, with each key value pair of the form (make sure to include spaces around the units):
        ::
            "<pipe diameter or area><' mm '/' m '/' mm^2 '/' m^2 '><any text you like>": <length of member in metres>, ...
    
    **Example:**
    ::
        gdn_sizedef = {
            "12.7 mm outer lab [m]": 3.10,
            "12.7 mm inner lab [m]": 6.55,
            "6.35 mm lab setup [m]": 1.40
        }
    """

    totalvol_m3 = 0

    for _, (key, value) in enumerate(pipelengths_dict.items()):

        if " mm " in key:
            pipediameter_m = float(key.split("mm")[0])/1e3
            pipelength_m = value
            totalvol_m3 += (np.pi * (pipediameter_m / 2) ** 2) * pipelength_m
        
        elif " m " in key:
            pipediameter_m = float(key.split("m")[0])
            pipelength_m = value
            totalvol_m3 += (np.pi * (pipediameter_m / 2) ** 2) * pipelength_m

        elif " mm^2 " in key:
            pipearea_m2 = float(key.split("mm^2")[0])/1e6
            pipelength_m = value
            totalvol_m3 += pipearea_m2 * pipelength_m
        
        elif " m^2 " in key:
            pipearea_m2 = float(key.split("m^2")[0])
            pipelength_m = value
            totalvol_m3 += pipearea_m2 * pipelength_m
    
    return totalvol_m3


def calc_partialpressures(ptarget_pa, qtarget, mdt_m3, gdn_m3, prefillpressure_pa=None):
    """Use this function to calculate the system pressure readings to accurately fill a gas mixture.

    **Parameters:**

    ptarget_pa
        float, the target initial pressure of the detonation (in Pascal).

    qtarget
        string, the composition by fraction of the gas mixture considered.

    mdt_m3
        float, the volume of the modular detonation tube.
    
    gdn_m3
        float, the approximate volume of a direct connection to the modular detonation tube.
    
    prefillpressure_pa
        float, the observed pressure of the system after settling from evacuation phase (in Pascal).
        Optional, defaults to None.

    **Example:**
    ::
        # Sizing parameters of the gas delivery network between gas bottles and MDT
        gdn_sizedef = {
            "12.7 mm outer lab [m]": 3.10,
            "12.7 mm inner lab [m]": 6.55,
            "6.35 mm lab setup [m]": 1.40
        }

        # Sizing parameters of the MDT
        mdt_sizedef = {
            "6200 mm^2 c-section [m]": 4.0
        }

        # Calculate volumes of the network
        gdnvolume_m3 = calc_networkvolume(gdn_sizedef)
        mdtvolume_m3 = calc_networkvolume(mdt_sizedef)

        # Calculate manometer readings an operator should see [Pa] for the piping system, when filling each gas
        _, sys_manometer_vals = calc_partialpressures(ptarget_pa=1e5, qtarget="C2H4:1 O2:3 N2:4", mdt_m3=mdtvolume_m3, gdn_m3=gdnvolume_m3, prefillpressure_pa=1000)
        print(sys_manometer_vals)

    Output:
    ::
        {'xx': 1000, 'C2H4': 12892.536570013475, 'O2': 48570.1462800539, 'N2': 100000.0}
    """

    labtemp_k = 298.15

    if prefillpressure_pa is None:
        # Imperfect vacuum pressure in tube and 13% of network |+| gas trapped in flexible tubing (can't evacuate 87% of the network)
        impurity_tube_pa = (800 * (0.13 * gdn_m3 + mdt_m3)  + 2e5 * (0.87 * gdn_m3)) / mdt_m3
        impurity_sys_pa = impurity_tube_pa * mdt_m3 / (mdt_m3 + gdn_m3)
    else:
        impurity_tube_pa = prefillpressure_pa * (mdt_m3 + gdn_m3) / mdt_m3
        impurity_sys_pa = prefillpressure_pa

    # Target moles in the tube, as given by ideal gas law, n = PV/RT
    targetmdtmoles = ptarget_pa * mdt_m3 / 8.314 / labtemp_k

    # Record cumulative totals
    qpp_tube_pa = {"xx": impurity_tube_pa}
    qrunningpp_sys_pa = {"xx": impurity_sys_pa}

    # Remove any gases from the component list if zero moles are used, and format q as dict
    qlist = [q.split(":") for q in qtarget.split(" ") if ":0" not in q]
    qdict = dict(zip([species[0] for species in qlist], [float(species[1]) for species in qlist]))
    qtotmoles = sum(list(qdict.values()))

    # Keep tabs of the cumulative pressure
    runningpp_tube_pa = impurity_tube_pa
    runningpp_sys_pa = impurity_sys_pa

    # Find the partial pressures for each of the species in the target gas composition
    for i, (species, qmoles) in enumerate(qdict.items()):

        # If the species is not last in filling, it uses extra volume from the gas network
        if i != len(qdict) - 1:
            mdtvolfraction = mdt_m3 / (mdt_m3 + gdn_m3)
            tube_deltapp_pa = (qmoles / qtotmoles) * ptarget_pa
            system_deltapp_pa = tube_deltapp_pa * mdtvolfraction
        # If the species is last in filling, it has to make up to the target pressure
        else:
            tube_deltapp_pa = ptarget_pa - runningpp_tube_pa
            system_deltapp_pa = ptarget_pa - runningpp_sys_pa
        
        # Record tube partial pressures
        runningpp_tube_pa += tube_deltapp_pa
        qpp_tube_pa[species] = tube_deltapp_pa

        # Record system partial pressures
        runningpp_sys_pa += system_deltapp_pa
        qrunningpp_sys_pa[species] = runningpp_sys_pa

    return qpp_tube_pa, qrunningpp_sys_pa


if __name__ == "__main__":

    # Sizing parameters of the gas delivery network between gas bottles and MDT
    gdn_sizedef = {
        "12.7 mm outer lab [m]": 3.10,
        "12.7 mm inner lab [m]": 6.55,
        "6.35 mm lab setup [m]": 1.40
    }

    # Sizing parameters of the MDT
    mdt_sizedef = {
        "6200 mm^2 c-section [m]": 4.0
    }

    # Calculate volumes of the network
    gdnvolume_m3 = calc_networkvolume(gdn_sizedef)
    mdtvolume_m3 = calc_networkvolume(mdt_sizedef)

    # Calculate manometer readings an operator should see [Pa] for the piping system, when filling each gas
    _, sys_manometer_vals = calc_partialpressures(ptarget_pa=1e5, qtarget="C2H4:1 O2:3 N2:4", mdt_m3=mdtvolume_m3, gdn_m3=gdnvolume_m3, prefillpressure_pa=800)
    print(sys_manometer_vals)

    pass
