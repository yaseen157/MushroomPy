<img src="./Documentation/img_src/Logo_1024x.png" alt="MushroomPy Logo" width=100>


# **Mushroom Py**

A collection of Python tools used in the processing of data from the 
[Virtual Test Facility](http://www.vtf.website/asc/wiki/bin/view/). 
This code is the product of a final-year student group project at the 
University of Southampton, to produce a modular rig for investigating 
detonation-combustion propagation behaviour.

*Authors: Yaseen Reza, Bruno del Mazo Canaleta*

## **Installation Pre-requisites**

> ### Python

MushroomPy is written for Python 3, and tested in Python 3.8 on 
Windows 10 (64-bit) machines. You should be able to check what 
version of Python you have by typing the following at the command
prompt

    $ py --version

If the prompt responds with `Python 3.x.x`, you have verified your 
installation - otherwise, it is recommended you download the official 
Python 3.8 distribution from [their website](https://www.python.org/downloads/).

> ### PIP

MushroomPy depends on a series of python modules, which are easily 
installed using the package installer for Python (pip). To check if 
you have pip, type the following at a command prompt

    $ pip --help

If the prompt responds, then you have verified your installation of
pip. Otherwise, follow the instructions [here](https://pip.pypa.io/en/stable/installing/). 
Note that if pip continues to be unresponsive following installation, 
it may be necessary to edit your *System Environment Variables* to 
include the directory that the pip executable is located in (this allows 
you to call pip from the command prompt, outside of its installed directory).
If you are missing any Python modules, they can be installed simply with 

    $ pip install module_name

> ### Shock and Detonation Toolbox

Class methods in the MushroomPy code frequently require the use of the 
[Shock and Detonation Toolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/).
It is required that you download and install both it and the [Cantera Software Package](https://cantera.org/) 
that it depends upon. Follow the [installation instructions](https://shepherd.caltech.edu/EDL/PublicResources/sdt/SDToolbox/sdt-install.pdf) 
given for Python installs by the toolbox webpage, which will walk you 
through the necessary steps.

To verify the toolbox and Cantera were installed correctly, begin
command line Python by typing the following at a command prompt

    $ py

Next, try importing the sdtoolbox and evaluating the following

    >>> import sdtoolbox
    >>> print(sdtoolbox.postshock.CJspeed(P1=100000, T1=298, q='C2H4:1 O2:3', mech='gri30_highT.cti'))

which if successful, will return the following Chapman-Jouget speed

    2373.202209779105


