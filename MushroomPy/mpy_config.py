# -*- coding: utf-8 -*-

"""
MushroomPy Module

Functions for interacting with the MushroomPy configuration file.

Updated January 2021
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
"""

__author__ = "Yaseen Reza"


import os
import tkinter as tk
import tkinter.filedialog      


def update_configfile(configdict):

    # Find the path of the MushroomPy configuration file
    module_directory = os.path.dirname(__file__)
    config_path = os.path.join(module_directory, "config.ini")

    # ADD LOGIC FOR DEFAULT CONFIGURATIONS HERE

    # dataroot: where the data for input/output operations will be kept
    # Check for key
    if "dataroot" not in configdict.keys():
        configdict["dataroot"] = ""
    # Check for value
    while os.path.exists(configdict["dataroot"]) is False:        
        # Wait for user acknowledgement of prompt
        input("Please specify a valid root directory in which MushroomPy produces and/or "
        "looks for data in - Press [Enter] to browse directories >")
        # Make a top level tkinter window and hide it
        root = tk.Tk()
        root.withdraw()
        # Bring hidden window to "focus"
        root.wm_attributes('-topmost', 1)
        # Populate the key
        configdict["dataroot"] = tkinter.filedialog.askdirectory(
            title="MushroomPy Setup - Choose Data Directory",
            parent=root,
            initialdir=os.path.expanduser("C:"))
        # Get rid of top level tkinter window
        root.destroy()
    
    # next configuration...
    # Check for key
    # ...
    # Check for value
    # ...

    # END OF LOGIC FOR DEFAULT CONFIGURATIONS

    # Configuration dictionary items are copied into the config file
    with open(config_path, "w+") as f:
        f.write("; MushroomPy Configuration\n")
        for _, (k, v) in enumerate(configdict.items()):
            f.write(f"{str(k)}={str(v)}\n")

    return

def parse_configfile():

    # Find the path of the MushroomPy configuration file
    module_directory = os.path.dirname(__file__)
    config_path = os.path.join(module_directory, "config.ini")

    # If the config file does not exist, create it
    if os.path.isfile(config_path) is False:
        with open(config_path, "w+") as f:
            pass
    
    # Parse the config file into a dictionary
    configdict = {}
    with open(config_path, "r") as f:
        for line in f.readlines():
            if "=" in line:
                set_k, set_v = line[:-1].split("=")
                configdict[set_k] = set_v

    # Update config.ini
    update_configfile(configdict=configdict)

    return configdict

if __name__ == "__main__":
    
    # Automatically perform config.ini setup
    parse_configfile()
