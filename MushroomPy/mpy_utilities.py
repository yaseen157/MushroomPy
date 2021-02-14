# -*- coding: utf-8 -*-

"""
MushroomPy Module

Common utilities for the MushroomPy package are contained within this module.

Updated January 2021
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
"""

__author__ = "Yaseen Reza"


import datetime


def gettimestr():
    """Call to return a string with the current time."""
    timenow = datetime.datetime.now()
    hour = timenow.hour
    minute = timenow.minute
    second = timenow.second
    ms = int(timenow.microsecond / 1000)

    return f"[{hour:02}:{minute:02}:{second:02}.{ms:03}]"


def getdatetimestr():
    """Call to return a string with the current date-time."""
    timenow = datetime.datetime.now()
    hour = timenow.hour
    minute = timenow.minute
    second = timenow.second
    date = timenow.date()
    return f"{str(date)}_{hour:02}-{minute:02}-{second:02}"

# Convert percentage dilution into molar fractions
def perc2mol(percent, fuelmols=1, oxmols=3):
    dilmoles = percent * (fuelmols+oxmols) / (100 - percent)
    return dilmoles


def mol2perc():
    return
