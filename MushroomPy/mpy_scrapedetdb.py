# -*- coding: utf-8 -*-

"""
MushroomPy Module

Web scraper for parsing data from the [detonation database](https://shepherd.caltech.edu/detn_db/html/data_sets.html).

Updated January 2021
Tested with:
    Python 3.8, Windows 10

Author Email: yr3g17@soton.ac.uk
"""

__author__ = "Yaseen Reza"


import requests
import pandas as pd
from bs4 import BeautifulSoup


ROOT_URL = "https://shepherd.caltech.edu/detn_db/"


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _getCellWidthURLext(queryfuel):
    """Given a queryfuel, return the URL extension of the ROOT_URL that provides Cell Width data."""

    # Get the datasets page and create a soup object of it
    page_datasets = requests.get(ROOT_URL+"html/"+"data_sets.html")
    soup = BeautifulSoup(page_datasets.content, "html.parser")
    # Search for "Cell Size" and keep the URL extension
    cellsize_urlelem = soup.find("a", string="Cell Size")
    # Make sure a valid element is detected before proceeding
    if cellsize_urlelem is None:
        return
    cellsize_urlext = cellsize_urlelem["href"]

    # Get the Cell size page and create a soup object of it
    page_cellsize = requests.get(ROOT_URL+"html/"+cellsize_urlext)
    soup = BeautifulSoup(page_cellsize.content, "html.parser")
    # Search for "Cell Width"+queryfuel and keep the URL extension
    cellwidth_urlelem = soup.find("a", string=f"Cell Width - {queryfuel} Fuel")
    # Make sure a valid element is detected before proceeding
    if cellwidth_urlelem is None:
        return
    cellwidth_urlext = cellwidth_urlelem["href"]

    return cellwidth_urlext


def _pandascleaner(pandasdf):
    """Given a dataframe scraped from the detonation database, return it with better formatting."""

    # The headers need to be reformatted
    headers = list(pandasdf.columns)
    for i in range(len(headers)):
        # If an undesirable character is detected in the header, left-strip and refactor the header
        while((headers[i][0] == "#") or (headers[i][0] == " ")):
            headers[i] = headers[i][1:]
        if "%" in headers[i]:
            headers[i] = "Percent".join(headers[i].split("%"))
        
        # The header should be properly case matched
        if "(" in headers[i]:
            splitheader = headers[i].split("(")
            headers[i] = "(".join([splitheader[0].title(), splitheader[1]])
        else:
            headers[i] = headers[i].title()

    # Refactor the dataframe column names with the reformatted versions
    pandasdf.columns = headers

    return pandasdf

def _recastaspandas(basepandadf, dicttorecast):
    """Given a dataframe and dictionary of initial conditions scraped from the database, return a full dataframe."""
    
    oldheaders_list = list(basepandadf.columns)
    newheaders_list = list(dicttorecast.keys())+oldheaders_list

    # Populate the new columns
    data_array = []
    for i_row in range(basepandadf.shape[0]):
        data_row = []
        for header in newheaders_list:
            # If the data comes from the dataframe
            if header in oldheaders_list:
                data_row.append(basepandadf[header][i_row])
            # Else the data came from the dictionary
            else:
                data_row.append(dicttorecast[header])
        data_array.append(data_row)

    # Construct a new dataframe
    returndf = pd.DataFrame(data_array, columns=newheaders_list)

    return returndf


def getCellWidthData(queryfuel):
    """Given a queryfuel, return a list of dictionaries of all Cell Width data scraped from the 
    detonation database.
    
    **Parameters:**
    queryfuel
        string, The text with which the website denotes a fuel you'd like to query.
    
    **Returns:**
    returndict
        dictionary, structured with a filename as a key, and its associated dataframe as a value.
    
    **Example:**
    ::
        c2h4data = getCellWidthData(queryfuel="C2H4")

        for _, (k, v) in enumerate(c2h4data.items()):
            print(f"This is the filename: | {k} | and here is the data:")
            print(v)
            print("")

    Output:
    ::
        ...

        This is the filename: | at172d.txt | and here is the data:
            Category  Fuel Sub-Category Oxidizer  Initial Pressure (kPa) Diluent Equivalence Ratio  Initial Temperature (K)  Cell Width (mm)
        0  cell size  C2H4        width      Air                   101.3                         1                   298.15             19.5
        1  cell size  C2H4        width      Air                   101.3                         1                   373.15             16.0
    
    **Example:**
    ::
        c2h4data = getCellWidthData(queryfuel="C2H4")
        oxidiserdict = {}

        for _, (k, v) in enumerate(c2h4data.items()):

            oxidisers = list(v["Oxidizer"])
            for oxidiser in oxidisers:
                if oxidiser in list(oxidiserdict.keys()):
                    oxidiserdict[oxidiser] += 1
                else:
                    oxidiserdict[oxidiser] = 0

        print("{'Oxidiser': Number of Hits} ==> ", oxidiserdict)
    
    Output:
    ::
        {'Oxidiser': Number of Hits} ==>  {'O2': 92, 'Air': 53}
    
    """

    urlext = _getCellWidthURLext(queryfuel=queryfuel)
    # Get the Cell width page and create a soup object of it
    page_cellwidths = requests.get(ROOT_URL+"html/"+urlext)
    soup = BeautifulSoup(page_cellwidths.content, "html.parser")
    # Create list of all available text files
    txtfiles = soup.find_all(string=lambda text: ".txt" in text)

    # Turn the list of discovered textfiles into a list of URL extensions
    urlext_list = []
    for txtfile in txtfiles:
        urlext = list(soup.find("a", string=txtfile)["href"])
        while (urlext.pop(0) != "/"):
            pass        
        urlext_list.append("".join(urlext))

    # Find all the relevant tables in the webpage and clean them up
    cellwidthdata_list = []
    for table in soup.find_all("table"):
        # The tables we are interested in have no attributes
        if len(table.attrs) == 0:

            # Keep all the detected table elements in a list
            tablestrings = []
            for string in table.strings:
                tablestrings.append(str(string))

            # Clean up the table elements
            tablestrings = list(filter(lambda element: element != "\n", tablestrings))
            for i in range(len(tablestrings)):
                # Remove leading '\n' characters
                while("\n" in tablestrings[i]):
                    tablestrings[i] = tablestrings[i][1:]
                # Remove leading and trailing spaces
                if len(tablestrings[i]) >= 2:
                    tablestrings[i] = (tablestrings[i][1:] if tablestrings[i][0] == " " else tablestrings[i])
                    tablestrings[i] = (tablestrings[i][:-1] if tablestrings[i][-1] == " " else tablestrings[i])
                # Remove trailing colons
                tablestrings[i] = (tablestrings[i][:-1] if tablestrings[i][-1] == ":" else tablestrings[i])
                # Cleaning up values
                if i%2 == 1:
                    # Check if value had units, and move units to the key
                    stringcomponents = tablestrings[i].split(" ")
                    if is_number(stringcomponents[0]) is True:
                        if stringcomponents[-1].isalpha() is True:
                            tablestrings[i-1] += f" ({stringcomponents[-1]})"
                            tablestrings[i] = float(stringcomponents[0])
                    # Check if value is a range, and mark for omission if it is
                    val= "".join(stringcomponents[:-1]) if stringcomponents[-1].isalpha() == True else "".join(stringcomponents)
                    if "-" in val:
                        if all(is_number(elem) for elem in val.split("-")) == True:
                            tablestrings[i-1] = "*delete"
                            tablestrings[i] = "*delete"

            # Delete keys and values marked for deletion
            tablestrings = list(filter(lambda element: element != "*delete", tablestrings))

            # Package every element of tablestrings in alternating dictionary keys and values
            cellwidthdata_list.append(dict(zip(tablestrings[::2], tablestrings[1::2])))

    # Package each dictionary of table data collected with the csv data that goes with it
    for i in range(len(cellwidthdata_list)):
        pandasdf = _pandascleaner(pd.read_table(f"{ROOT_URL}{urlext_list[i]}", delimiter=",", dtype=float))
        cellwidthdata_list[i] = _recastaspandas(basepandadf=pandasdf, dicttorecast=cellwidthdata_list[i])
    returndict = dict(zip(txtfiles, cellwidthdata_list))

    return returndict


if __name__ == "__main__":

    pass
