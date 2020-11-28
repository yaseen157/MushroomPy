__author__ = "Yaseen Reza"

import copy
import datetime
import math
import os
import re
import statistics as st
import time
import warnings
from itertools import cycle

# import openpyxl (this or a similar excel writing engine must be installed to export XLSX reports)
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sdtoolbox.postshock import CJspeed, PostShock_fr, PostShock_eq


# Email: yr3g17@soton.ac.uk


def _gettimestr():
    """Global function, call to return a string with the current time."""
    timenow = datetime.datetime.now()
    hour = timenow.hour
    minute = timenow.minute
    second = timenow.second
    ms = int(timenow.microsecond / 1000)

    return f"[{hour:02}:{minute:02}:{second:02}.{ms:03}]"


def _getdatetimestr():
    """Global function, call to return a string with the current date-time."""
    timenow = datetime.datetime.now()
    hour = timenow.hour
    minute = timenow.minute
    second = timenow.second
    date = timenow.date()
    return f"{str(date)}_{hour:02}-{minute:02}-{second:02}"


class SDT0d:
    """**Zero-dimensional** detonation analysis tool, using methods from the **Shock & Detonation Toolbox**."""

    def __init__(self, p_initial, t_initial, q_initial):
        self.P1 = p_initial
        self.T1 = t_initial
        self.q = q_initial
        self.mech = 'gri30_highT.cti'

    def calc_cjspeed(self):
        # Calculation for the Chapman-Jouget speeds
        cj_speed_mps = CJspeed(P1=self.P1, T1=self.T1, q=self.q, mech=self.mech, fullOutput=False)

        return cj_speed_mps

    def calc_cjpressure(self):
        # Calculation for the Chapman-Jouget speeds
        cj_speed_mps = self.calc_cjspeed()

        # Post shock (equilibrium) pressure to be determined
        gascj = PostShock_eq(U1=cj_speed_mps, P1=self.P1, T1=self.T1, q=self.q, mech=self.mech)
        cj_pressure_pa = gascj.P

        return cj_pressure_pa

    def calc_vnpressure(self):
        # Calculation for the Chapman-Jouget speeds
        cj_speed_mps = self.calc_cjspeed()

        # Post shock (frozen) pressure to be determined
        gasvn = PostShock_fr(U1=cj_speed_mps, P1=self.P1, T1=self.T1, q=self.q, mech=self.mech)
        vn_pressure_pa = gasvn.P

        return vn_pressure_pa


class AMROC1d:
    """**One-dimensional** detonation analysis tool, intended for use with simulation data from
    '../amroc/clawpack/applications/euler_chem/1d/DetonationTube'. This class is designed for use with a
    fuel:oxidiser:monodiluent setup.

    Where the root directory contains this module (root/mdsimtools.py), simulation text files should be stored here:
        root/1D_AMROC_data/<driver_kpa>_<driven_kpa>_<moldiluent-unitymolesfuel_ratio>/<diluent>_<mechanism>/dat_xxx.txt
    For example, 'root/1D_AMROC_data/50_10_2_data/Ar_C2H4_Jachi/dat_134.txt'.

     """

    def __init__(self):
        # Where is this Python tools file located?
        # Python module directory
        self.module_path = os.path.dirname(__file__)
        # Python module name
        self.module = os.path.basename(__file__).replace(".py", "")
        # Python class name
        self.classname = type(self).__name__

        # Where are the Input and Output directories?
        # Output data directory (e.g. ...\_output\hpcsimtools\AMROC1d\...)
        self.output_path = os.path.join(self.module_path, "_output", self.module, self.classname)
        # Input data directory
        self.simdata_path = os.path.join(self.module_path, "1D_AMROC_data")

        # Check if the simulation data folder is available, else make it
        if not os.path.exists(self.simdata_path):
            warnmsg = f"Simulation data path does not exist, generating a new one ({self.simdata_path})"
            warnings.warn(warnmsg)
            os.makedirs(self.simdata_path)

        # Find all available data folders in the simulation data folder
        datadirlist = [datafile for datafile in os.listdir(self.simdata_path) if "data" in datafile]

        # Produce an ordered list of data case folders, if any were discovered
        datadir_sorted = self._getsortedcataloguekeys(cataloguekeyslist=datadirlist)

        # Construct a file hierarchy tree/catalogue
        datacatalogue = dict(zip(datadir_sorted, [None for _ in range(len(datadir_sorted))]))

        # For each data folder, add tree branches for mechanism type
        for casename in datadir_sorted:
            case_path = os.path.join(self.simdata_path, casename)
            mechanisms = os.listdir(case_path)
            datacatalogue[casename] = dict(zip(mechanisms, [None for _ in range(len(mechanisms))]))

            # For each mechanism type, add tree branches for data files available
            for _, (mechanism, _) in enumerate(datacatalogue[casename].items()):
                mechanism_path = os.path.join(self.simdata_path, casename, mechanism)
                datafiles = [file for file in os.listdir(mechanism_path) if "dat_" in file]

                # Rearrange the file list sorted numerically rather than alphabetically
                filenumbers = sorted([int(re.findall("\d+", element)[0]) for element in datafiles])
                datafiles_sorted = [f"dat_{element}.txt" for element in filenumbers]
                datacatalogue[casename][mechanism] = datafiles_sorted

        # A list of every single data folder, the mechanism/diluent combination detected, and the associated data files
        self.datacatalogue = datacatalogue

        # As we read lines of data, store it in a volatile catalogue as required
        self.volatilecatalogue = {}

        # As statistics become available, store them in here
        self.statistics = {}

    def case_read(self, datacasefolder):

        # Is the folder name valid? Proceed if the folder can be found in the catalogue
        if datacasefolder in self.datacatalogue.keys():

            # Inform the user task is beginning
            print(f"{_gettimestr()} Reading case '{datacasefolder}'...")

            # Create an entry in the volatile data catalogue
            self.volatilecatalogue[datacasefolder] = {}

            # For each set of data contained by a mechanism folder in the data case folder (described by full catalogue)
            for _, (mechanism, textfiles) in enumerate(self.datacatalogue[datacasefolder].items()):

                temptxtdict = {}

                # For each text file in the list of text files within a mechanism
                for textfile in textfiles:
                    datatxt_path = os.path.join(self.simdata_path, datacasefolder, mechanism, textfile)

                    # Produce a list of headers from the datafile
                    header = pd.read_csv(datatxt_path, nrows=0, sep=',\s+', skipinitialspace=True, engine='python')
                    header.columns = header.columns.str.replace(' ', '')
                    headers = header.columns[0].split(";")

                    # Produce a dataframe of all the data read
                    data = pd.read_csv(datatxt_path, skiprows=1, delimiter=" ", names=headers[1:]).astype(float)
                    data.drop(data.columns[-1], axis=1, inplace=True)

                    # Package nicely into pandas dataframes in the master data dictionary, and sort by timestamp
                    temptxtdict[textfile] = {"Time Elapsed": float(headers[0]), "DataframeObj": data}
                    temptxtdict = (
                        {k: v for k, v in
                         sorted(temptxtdict.items(), key=lambda kvtuple: kvtuple[1]["Time Elapsed"])})
                    self.volatilecatalogue[datacasefolder][mechanism] = temptxtdict

        # The folder name queried could not be found in the available data catalogue, return a warning
        else:
            warnmsg = f"Could not find '{datacasefolder}' in available data cases!"
            warnings.warn(warnmsg)

    def case_remove(self, datacasefolder):

        # Inform the user task is beginning
        print(f"{_gettimestr()} Removing data case '{datacasefolder}'...")

        # Is the folder name valid? Proceed if the folder can be found in the catalogue
        if datacasefolder in self.volatilecatalogue.keys():
            removed = self.volatilecatalogue.pop(datacasefolder)
            return removed
        # The folder name queried could not be found in the volatile data catalogue, return a warning
        else:
            warnmsg = f"Could not find '{datacasefolder}' in loaded data cases!"
            warnings.warn(warnmsg)

    def data_headerscatalogue(self):

        # Preamble
        volcat = self.volatilecatalogue

        # Build a dictionary to catalogue available headers
        headers_dict = {}

        # For each datacase available in the volatile catalogue
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
            mechanisms = volcat[datacase].keys()

            # Initialise a dictionary to store a catalogue of volatile headers dict={datafolder : {mechanism: {}}}
            headers_dict[datacase] = dict(zip(mechanisms, [{} for _ in range(len(mechanisms))]))

            # For each mechanism in a data case, catalogue the available headers to a list
            for mechanism in mechanisms:
                headers_dict[datacase][mechanism] = list(
                    self.volatilecatalogue[datacase][mechanism]["dat_0.txt"]["DataframeObj"].columns)

        return headers_dict

    def data_statistics(self):

        # Preamble, before user is notified task is running
        volcat = self.volatilecatalogue
        statistics_dict = self.statistics

        # Check if statistics have already been generated, and return if the statistics dictionary is already up to date
        missingcases_set = volcat.keys() - statistics_dict.keys()
        if len(missingcases_set) == 0:
            return self.statistics

        # Further preamble, before user is notified task is running
        headers_dict = self.data_headerscatalogue()

        # Inform the user task is beginning and initialise variables to track the progress metric of this method
        print(f"{_gettimestr()} Generating statistics...")
        start = time.perf_counter()
        casecount_current = 0
        casecount_target = str(len(missingcases_set))

        # For each datacase available in the volatile catalogue
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
            mechanisms = volcat[datacase].keys()

            # The datacase only needs to be read if its missing from the statistics catalogue
            if datacase in missingcases_set:

                # Inform the user task is progressing
                casecount_current += 1
                casecount_str = format(casecount_current, "0" + str(len(casecount_target))) + "/" + casecount_target
                print(f"{_gettimestr()} ...for case {casecount_str} | '{datacase}'")

                # Flesh out the statistics dictionary with more branches
                statistics_dict[datacase] = dict(zip(mechanisms, [{} for _ in range(len(mechanisms))]))

                # For each mechanism in a data case, find the available headers that aren't position related
                for mechanism in mechanisms:
                    headers = [header for header in headers_dict[datacase][mechanism] if header != "x"]

                    # Initialise parameters to build a pandas dataframe to track the combustion front
                    combustiontracker = [["time (s)", "x"]]
                    _, _, q, _ = self.case_initialconditions(datacase=datacase, mechanism=mechanism)
                    flamefront_idx = 0

                    # For each header
                    for header in headers:

                        # Build a pandas data frame with these headers
                        rows = [["time (s)", "maximum", "mean", "minimum"]]

                        # For every time step (This takes the longest)
                        for txtfile in volcat[datacase][mechanism].keys():
                            txtfiledata = volcat[datacase][mechanism][txtfile]
                            # Read data from the volatile catalogue
                            timestamp, df = txtfiledata.values()
                            # Gather all the data from a header, and return statistical metadata
                            y = [element for element in list(df[header])]
                            rows.append([timestamp, max(y), sum(y) / len(y), min(y)])

                            # If the header contains a species name and matches the first species described in q
                            if header.split("-")[-1] == q.split(":")[0]:
                                # Find and record the position of the reaction front against time
                                try:
                                    flamefront_idx = y.index(rows[1][1])
                                    combustiontracker.append([timestamp, df["x"][flamefront_idx]])
                                except ValueError:
                                    combustiontracker.append([timestamp, max(df["x"])])

                        # Generate and save the pandas dataframes for generic headers
                        tempdf = pd.DataFrame(rows[1:], columns=rows[0])
                        statistics_dict[datacase][mechanism][header] = tempdf

                    # Generate and save the pandas dataframes for the combustion front
                    tempdf = pd.DataFrame(combustiontracker[1:], columns=combustiontracker[0])
                    statistics_dict[datacase][mechanism]["ReactionFront"] = tempdf

        # Update the class' statistics definition
        self.statistics = statistics_dict

        # End a timer and return the time taken for the statistics method to run
        finish = time.perf_counter()
        print(f"{_gettimestr()} (Generated statistics in {(finish - start):.1f} s)")

        return self.statistics

    def export_statistics(self):

        # Preamble, before user is notified task is running
        volcat = self.volatilecatalogue
        headers_dict = self.data_headerscatalogue()
        statistics_dict = self.data_statistics()

        # Inform the user task is beginning
        print(f"{_gettimestr()} Exporting statistics...")

        # If the output directory being exported doesn't exist, make it
        outputdir_path = os.path.join(self.output_path, "statsCSVs")
        if not os.path.exists(outputdir_path):
            os.makedirs(outputdir_path)

        # Track the progress of the function using a file-counting variable
        filecount = 0

        # For every data case, every mechanism within, and every header there-in, produce CSVs
        # For each datacase available in the volatile catalogue
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
            mechanisms = volcat[datacase].keys()

            # For each mechanism in a data case, find the available headers that aren't position related
            for mechanism in mechanisms:
                headers = [header for header in headers_dict[datacase][mechanism] if header != "x"]

                # For each header
                for header in headers:
                    outputfile_name = f"{datacase}+{mechanism}+{header}.csv"
                    outputfile_path = os.path.join(outputdir_path, outputfile_name)

                    statistics_dict[datacase][mechanism][header].to_csv(outputfile_path, index=False, header=True)
                    filecount += 1

        # Return the total number of counted files exported
        print(f"{_gettimestr()} >> Exported {filecount} file(s) to '{outputdir_path}'.")

        return

    def plot_mechanisms(self, colourseed=None):

        # Preamble, before user is notified task is running
        volcat = self.volatilecatalogue
        statistics_dict = self.data_statistics()
        simscores_dict = self._getsimilarityscores()

        # Inform the user task has proceeded
        print(f"{_gettimestr()} Plotting mechanism comparisons...")

        # If the output directory being exported doesn't exist, make it
        outputdir_path = os.path.join(self.output_path, "statsPNGs")
        if not os.path.exists(outputdir_path):
            os.makedirs(outputdir_path)

        # Track the progress of the function using a file-counting variable
        filecount = 0

        # The list of headers to plot are headers that are common across mechanisms, not relating to species or position
        headers_set = self._getcommonheaders()
        headers = [header for header in headers_set]

        # For every data case
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
            mechanisms = volcat[datacase].keys()

            # For the data case considered, determine the set of diluents made available in the mechanisms
            mechanism_parts = [dil_mechanism.split("_") for dil_mechanism in mechanisms]
            diluents = set([sublist[0] for sublist in mechanism_parts])

            # The plots can be coloured by mechanism, and can use a time dependent or independent colour seed
            cycled_clrs = cycle(list(['limegreen', 'darkorchid', 'indigo', 'darkgoldenrod',
                                      'royalblue', 'darkslategrey', 'orange', 'darkturquoise', 'forestgreen',
                                      'red', 'darkviolet', 'yellowgreen', 'plum', 'navy',
                                      'gold', 'dodgerblue', 'teal', 'darkred']))
            if colourseed:
                for _ in range(datetime.datetime.now().second % 7):
                    next(cycled_clrs)

            # Each mechanism should be mapped by a colour dictionary to a colour to ensure uniformity
            clr_dict = {}
            for _ in mechanisms:
                clr_dict = dict(zip(mechanisms, [next(cycled_clrs) for _ in range(len(mechanisms))]))

            # The names and line styles of the plot y-series
            ydata_keys = ["maximum", "mean"]
            ls_dict = dict(zip(ydata_keys, ["-", "--"]))

            # For each header
            for header in headers:

                # For each diluent, superimpose mechanism plots and save them
                for diluent in diluents:

                    # Determine the similarity scores for peak and average values, across mechanisms with common diluent
                    simscore1 = st.fmean(simscores_dict[datacase][header]["maximum"][diluent])
                    simscore2 = st.fmean(simscores_dict[datacase][header]["mean"][diluent])

                    # Figure preamble
                    fig = plt.figure(figsize=[8, 8], dpi=100)
                    ax = fig.add_axes([0.14, 0.12, 0.72, 0.8])  # left bottom width height
                    ax.grid(True)
                    ax.set_title(f'{header} (S1={simscore1:.4f}; S2={simscore2:.4f})')
                    ax.set_xlabel("time (s)")
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

                    # For each mechanism available in the data case, now try to create superimposed plots
                    for mechanism in mechanisms:
                        # Mechanisms are separated by diluent
                        if diluent in mechanism:
                            xdata = statistics_dict[datacase][mechanism][header]["time (s)"]
                            # Plot the individual series
                            for ydata_key in ydata_keys:
                                ydata = statistics_dict[datacase][mechanism][header][ydata_key]
                                label = f"{mechanism} ({ydata_key})"
                                ax.plot(xdata, ydata, label=label, c=clr_dict[mechanism], ls=ls_dict[ydata_key])

                    # Finish the plot with a legend
                    ax.legend()

                    # Save the resulting figure
                    outputfile_name = f"{diluent}+{datacase}+{header}.png"
                    outputfile_path = os.path.join(outputdir_path, outputfile_name)
                    plt.savefig(fname=outputfile_path)
                    plt.close('all')

                    # Increment the file-counting metric
                    filecount += 1

        # Return the total number of counted files exported
        print(f"{_gettimestr()} >> Exported {filecount} file(s) to '{outputdir_path}'.")

    def _getcommonheaders(self):
        """In order to compare several mechanisms, only compare headers that are common between them."""

        # Only proceed if there are data cases present in the volatile data catalogue
        if len(self.data_headerscatalogue()) > 0:

            # The list of headers to plot should be obtained from headers that are common across mechanisms
            headersets_list = []
            for _, (datacase, mechanisms) in enumerate(self.data_headerscatalogue().items()):

                # Take from index [1:] to omit 'x' (positional data) from common headers
                headersets_list = [set(headerlist[1:]) for headerlist in mechanisms.values()]

                # If a common header could not be found, it is likely no mechanisms exist, or mechanism data is empty
                if len(headersets_list) == 0:
                    raise ValueError(f"Could not find common header, check data integrity in '{datacase}'.")

            # The set of common headers can be found from set intersection of a nested list of volatile headers
            intersectingheaders_set = set.intersection(*headersets_list)

            # Even if they prove to be common headers, remove any headers relating to gaseous species molar fractions
            commonheaders_list = [elem for elem in intersectingheaders_set if "Y_" not in elem]

        # If no data is loaded in the volatile catalogue, then the list of available headers is simply empty
        else:
            commonheaders_list = []

        return commonheaders_list

    def _getcommontimestep(self):
        """In order to appropriately compare two mechanisms, find the common time steps between them."""

        # Store all available timestamps in a nested list, initialise the list here
        tslist_nested = []

        # Add all known timestamps to the nested list
        for _, (datacase, mechanisms) in enumerate(self.volatilecatalogue.items()):
            for mech in mechanisms:
                timesteps = [float(mechanisms[mech][txtfile]["Time Elapsed"]) for txtfile in mechanisms[mech].keys()]
                tslist_nested.append(set(timesteps))

        # A set intersection reveals the list of timesteps common amongst all datafiles queried
        intersectingtimesteps = sorted(list(set.intersection(*tslist_nested)))

        return intersectingtimesteps

    def _getsimilarityscores(self):
        """To compare several mechanisms for a datacase, a self-proposed method is used to compare the results of
        mechanisms, by assigning a similarity score (Si) to each mechanism.

        S1 = similarity score for the peak header values from a mechanism.
        S2 = similarity score for the mean header values from a mechanism.

        For the curve (x1, y1) of peak values produced by mechanism_1, to the curve (xn, yn) of peak values described
        by mechanism_n, we seek a performance metric that allows for the simple evaluation of mechanism agreement. The
        similarity comparison can only happen if the number of elements in x1 is equal to curves up to xn.

        A similarity score S for the curve (x1, y1) at a *particular* value of x (x') can be computed as follows:
        mu' = average(y'1 ... y'n)
        S' = 1 / exp(absolutevalue((y'1 / mu') - 1))

        The similarity score has the attractive properties that y-data is normalised by the average of the y-data and
        is therefore independent of scale-factors or unit conversions, moreover S' is bounded by [0, 1], allowing for
        dimensionless comparisons. Because this function scales as an exponential of Euler's number, a score is easily
        predicted if the average distance of a curve to the curve mean is known.
        """

        # Preamble, find common headers and timesteps with which a similarity score can be computed
        headers_isct = self._getcommonheaders()
        timesteps = self._getcommontimestep()
        statistics_dict = self.data_statistics()

        # Set up a list of template dictionaries for scores to be copied
        scorekeys = ["maximum", "mean"]
        # e.g. score_dicts = [{header : {"maximum" : {}, "mean" : {}}}, ..., {...}]
        score_dicts = [{header: dict(zip(scorekeys, [{} for _ in range(len(scorekeys))]))} for header in headers_isct]

        # Create a dictionary to store all similarity scores in
        # e.g. dict={datacase : {header : {"maximum" : {}, "mean" : {}}, header2 : {...}}, datacase2 : {...}, ...}
        similarityscore_dict = {}
        for datacase in statistics_dict.keys():
            similarityscore_dict[datacase] = {}
            for dictionary in score_dicts:
                similarityscore_dict[datacase].update(copy.deepcopy(dictionary))

        # Function for finding the similarity score [a, b, c] for nested lists [[a1, b1, c1][...][an, bn, cn]]
        def find_similarityscores(nestednumberslists):

            # The nested number lists need to have their 'nth' elements averaged into a new list
            np_data = np.array(nestednumberslists)
            np_mu = np.average(np_data, axis=0)

            # Calculate and append the similarity score to a results list of identical size to input list
            scores = []
            for numberlist in nestednumberslists:
                with np.errstate(divide='ignore', invalid='ignore'):
                    denominator = np.exp(np.absolute(np.divide(np.array(numberlist), np_mu) - 1))
                denominator[np.isnan(denominator)] = 1
                scores.append(list(np.divide(1, denominator)))

            # All the nested list inputs have been converted into similarity scores, take the average to make one list
            np_averagescores = np.average(np.array(scores), axis=0)

            return np_averagescores.tolist()

        # For each datacase as passed in with the statistics dictionary
        for _, (datacase, mechanisms) in enumerate(statistics_dict.items()):

            # Step 1: Extract and aggregate statistical data for each datacase
            # Construct a temporary dictionary for a datacase
            headerstats_dict = {}
            headerstats_dict.clear()
            for dictionary in score_dicts:
                headerstats_dict.update(copy.deepcopy(dictionary))

            # All the statistical data from all mechanisms should be packaged into the datacase-specific temp dictionary
            for mechanism in mechanisms:

                # Identify the diluent used in the mechanism - we need to separate scored curves by the diluent used
                diluent = mechanism.split("_")[0]

                # For every header, prepare its statistical data
                for header in headers_isct:
                    headerdf = statistics_dict[datacase][mechanism][header]

                    # For every header statistic (maximum, mean), create an empty list to append nested lists of data to
                    for key in scorekeys:
                        # The nested list parent only needs to be created if it doesn't already exist
                        if diluent not in list(headerstats_dict[header][key].keys()):
                            headerstats_dict[header][key][diluent] = []

                    # Record the average value for each mechanism-unique header max and mean
                    for key in scorekeys:
                        # Since different data sets may have different time steps, only append data from common steps
                        values = [float(headerdf[key][i]) for i in range(len(headerdf[key])) if
                                  float(list(headerdf["time (s)"])[i]) in timesteps]
                        # Time step synchronsised data can be stored in nested lists
                        headerstats_dict[header][key][diluent].append(values)

            # Step 2: Compute the similarity scores and store this data
            # Enter a loop hierarchy similar to that from Step 1
            for mechanism in mechanisms:
                diluent = mechanism.split("_")[0]
                for header in headers_isct:
                    for key in scorekeys:

                        # Try to determine the similarity score for time synchronised curves
                        try:
                            similarityscore_dict[datacase][header][key][diluent] = find_similarityscores(
                                headerstats_dict[header][key][diluent])
                        # If a TypeError was encountered, it's most likely because some elements from txt files read nan
                        except TypeError:
                            pathstring = f"...\\{datacase}\\{mechanism}"
                            raise TypeError(f"Could not find similarity score, check data integrity in '{pathstring}'.")

        return similarityscore_dict

    def data_detonationcheck(self):
        """Categorise the input data to determine the type of combustion taking place. There are two different checks
        available in this method.

        Check 1: pressure, what is the average pressure rise in the tube? Used to check simulation health.
        Check 2: flamefront, what is the speed of the combustion front? Differentiates between combustion types.
        """

        # Preamble
        statistics_dict = self.data_statistics()
        volcat = self.volatilecatalogue

        # Initialise a dictionary to store the results in
        detonation_dict = {}

        def detonationcheck_pressure(statistics_mechanism_df):
            # The initial pressure is the minimum pressure in the tube
            press_i = min(list(statistics_mechanism_df["Pressure"]["minimum"]))
            press_fbar = st.fmean(list(statistics_mechanism_df["Pressure"]["mean"]))

            # If the ratio of final tube mean pressure to initial pressure returns a nan, the simulation likely failed
            if np.isnan(press_fbar / press_i):
                return "Error"
            else:
                return None

        def detonationcheck_flamefront(statistics_mechanism_df):
            # What is the furthest distance along the tube reached (in centimetres)
            furthestdistance_x = max(statistics_mechanism_df["ReactionFront"]["x"])
            # What is the index this distance is first reached
            furthestdistance_idx = list(statistics_mechanism_df["ReactionFront"]["x"]).index(furthestdistance_x)
            # What is the time taken to reach this distance (in seconds)
            furthestdistance_s = statistics_mechanism_df["ReactionFront"]["time (s)"][furthestdistance_idx]

            # What is the index of the final item in the progression tracker of the flame front?
            listfinalitem_idx = len(statistics_mechanism_df["ReactionFront"]["x"]) - 1

            # If the furthest distance travelled index is equal to than the total tracked index, the sim did not finish
            if listfinalitem_idx == furthestdistance_idx:
                return "DNF"

            # Else, work out the flame front average Mach number
            else:
                # What is the average speed of the wave at this point?
                avgwavespeed_mps = furthestdistance_x / 100 / furthestdistance_s
                # What is the speed of sound?
                a_mps = st.fmean(list(statistics_mechanism_df["SpeedofSound"]["mean"])[0:furthestdistance_idx])

                mach = avgwavespeed_mps / a_mps
                if mach >= 1:
                    return True
                else:
                    return False

        # For each datacase available in the volatile catalogue
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
            mechanisms = volcat[datacase].keys()

            # Each datacase should have a dictionary to store the mechanism results in
            detonation_dict[datacase] = {}

            # For each mechanism in the mechanisms available to the datacase, find the mechanism's header dataframes
            for mechanism in mechanisms:
                mechdfs = statistics_dict[datacase][mechanism]

                # Check data for nan simulation errors, else label the data case and mechanism with detonation status
                if detonationcheck_pressure(statistics_mechanism_df=mechdfs) == "Error":
                    detonation_dict[datacase][mechanism] = "Error"
                else:
                    detonation_dict[datacase][mechanism] = detonationcheck_flamefront(statistics_mechanism_df=mechdfs)

        return detonation_dict

    def _getsortedcataloguekeys(self, cataloguekeyslist):

        # Find all available data folders in the simulation data folder
        datadirlist = cataloguekeyslist

        # Produce an ordered list of data case folders, if any were discovered
        datadir_sorted = []
        if len(datadirlist) > 0:
            # Break up all the data folder names into lists containing [driver_pa, driven_pa, diluent_molfrac, "data"]
            datadir_parts = sorted([re.split("_", element) for element in datadirlist])
            # Group driver/driven pressures, and diluent moles into sets (but of type list)
            datadir_nums = [sorted(list({int(sublist[x]) for sublist in datadir_parts})) for x in
                            range(len(datadir_parts[0]) - 1)]
            # Reconstruct the original folder names
            for driverpressure in datadir_nums[0]:
                for drivenpressure in datadir_nums[1]:
                    for diluentmoles in datadir_nums[2]:
                        testcase = [str(driverpressure), str(drivenpressure), str(diluentmoles)]
                        datadir_sorted += ["_".join(sublist) for sublist in datadir_parts if sublist[:-1] == testcase]

        return datadir_sorted

    def export_detonationreport(self):

        # The report generator is dependent on a module not required by the rest of the code
        try:
            import openpyxl
        except ImportError:
            raise ImportError("Could not create report, ensure module 'openpyxl' is installed!")

        # Preamble, find the cases of detonations
        detonationcheck_dict = self.data_detonationcheck()

        # Inform the user task is beginning
        print(f"{_gettimestr()} Exporting detonation booleans...")

        # If the output directory being exported doesn't exist, make it
        outputdir_path = os.path.join(self.output_path, "reports")
        if not os.path.exists(outputdir_path):
            os.makedirs(outputdir_path)

        # For all the mechanisms in the detonationcheck_dict, find the set of mechanisms present
        mechanisms = list({j for i in detonationcheck_dict.values() for j in i})

        # Create a new dictionary with the intention of storing detonation by the mechanism used
        mechanismdetonation_dict = dict(zip(mechanisms, [{} for _ in mechanisms]))

        # Define a function which will be used later by the export engine to colour code the output excel worksheet
        def highlight_bools(val):
            # Successful detonations
            if val is True:
                colour = '#7FFF00'  # CSS chartreuse
            elif val == "DNF":
                colour = '#00FA9A'  # CSS mediumspringgreen
            elif val == "Weak":
                colour = '#7FFFD4'  # CSS aquamarine
            # Unsuccessful detonations
            elif val is False:
                colour = '#FFD700'  # CSS darkorange
            # Data Error encountered
            elif val == "Error":
                colour = '#FF0000'  # CSS red
            elif val == "n/a":
                colour = '#E6E6FA'  # CSS lavender
            # Other cells
            else:
                return ""
            return f"background-color: {colour}"

        # Step 1: Refactor the detonation boolean dictionary by mechanism used
        # For each mechanism
        for mechanism in mechanisms:

            # Check if the mechanism was in the datacase, and if it was, append to the refactored dictionary
            for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(self.volatilecatalogue.keys())):
                if mechanism in detonationcheck_dict[datacase].keys():
                    mechanismdetonation_dict[mechanism][datacase] = detonationcheck_dict[datacase][mechanism]

        # Step 2: For each mechanism, export a file containing a boolean of all the detonation statuses
        outputfile_name = f"detonations_{_getdatetimestr()}.xlsx"
        outputfile_path = os.path.join(outputdir_path, outputfile_name)
        writer = pd.ExcelWriter(outputfile_path, engine='openpyxl')

        # For each mechanism
        for mechanism in mechanisms:

            # For each datacase contained within the mechanism
            datacases = list(set(mechanismdetonation_dict[mechanism].keys()))
            datacases_parts = sorted([datacasename.split("_") for datacasename in datacases], key=lambda x: int(x[0]))

            # Find the set of diluents that were at any point used with the mechanism being investigated
            diluentslist = sorted(list({int(sublist[2]) for sublist in datacases_parts}))

            # For data folder structure 'kpa_kpa_n_data', produce an ordered set of 'kpa_kpa' "roots"
            kpa_kpa = list({"_".join(sublist[0:2]) for sublist in datacases_parts})
            kpa_kpa = ["_".join(x) for x in
                       sorted([p.split("_") for p in kpa_kpa], key=lambda x: (int(x[0]), int(x[1])))]

            # Build a pandas data frame with these headers
            rows = [["driver_driven"] + diluentslist]

            # For the available pressure cases in the given mechanism add the detonation booleans to the dataframe data
            for pressurecase in kpa_kpa:

                # Initialise an empty data row to be appended to 'rows' later
                datarow = []

                # Produce a list of datacase names, not all of which are valid (missing data for missing diluent cases)
                for diluentmole in diluentslist:
                    datacase_tocheck = "_".join([pressurecase, str(diluentmole), "data"])

                    # If the data folder exists for the given diluent quantity, append the detonation check boolean
                    if datacase_tocheck in datacases:
                        datarow.append(mechanismdetonation_dict[mechanism][datacase_tocheck])
                    # If the data folder did not exist, append a "not applicable" empty value placeholder
                    else:
                        datarow.append("n/a")

                # Finally, the pressure case tested, along with all the diluent moles tested, is appended to 'rows'.
                rows.append([pressurecase] + datarow)

            # Generate the pandas dataframes and store into excel sheets
            df = pd.DataFrame(rows[1:], columns=rows[0])
            df.style.applymap(highlight_bools).to_excel(writer, index=False, engine="openpyxl", sheet_name=mechanism)

        # Save excel worksheet and write-out
        writer.save()
        print(f"{_gettimestr()} >> Exported Report '{outputfile_name}' to '{outputdir_path}'.")

    def case_initialconditions(self, datacase, mechanism):

        # Preamble
        volcat = self.volatilecatalogue

        # Work out the number of moles of diluent (stored in the name of the data case)
        diluentmoles = float(datacase.split("_")[2])

        # Set up a simple dataframe to refer back to, of initial conditions
        initialdf = volcat[datacase][mechanism]["dat_0.txt"]["DataframeObj"]

        # First find initial conditions P1, T1, and q (follows)
        p1 = st.mode(initialdf["Pressure"])
        t1 = st.mode(initialdf["Temperature"])

        # Create a dictionary of the species present in the mechanism, along with the species fraction
        species = sorted([y_0x.split("-") for y_0x in initialdf.columns if "Y_" in y_0x], key=lambda x: x[1])
        specieslist = ["-".join(sublist) for sublist in species]
        speciesfracs = [st.mode(initialdf[species]) for species in specieslist]
        speciesdict = dict(zip([species.split("-")[1] for species in specieslist], speciesfracs))

        # Work out the fraction of the mixture that is diluent, if at all
        diluentfrac = 1 - sum(speciesfracs)

        # If diluent is present, it is necessary to determine the proper molar proportions
        if diluentfrac > 0:
            summoles = diluentmoles / diluentfrac
        else:
            summoles = 1

        # Create a string to describe the gas mixture, that the Shock and Detonation Toolbox can also understand
        q = " ".join([f"{k}:{v * summoles:.2f}" for _, (k, v) in enumerate(speciesdict.items()) if v > 0])

        # If diluent is present, add it to the string
        if diluentfrac > 0:
            # Work out the diluent gas (stored in the name of the mechanism)
            diluent = mechanism.split("_")[0]
            q += f"{diluent}:{diluentmoles}"

        return p1, t1, q, diluentmoles

    def data_vn2cjratio(self, fast=True):

        # Preamble
        volcat = self.volatilecatalogue

        # Initialise a dictionary to store the results in
        znd_pratio_dict = {}

        # For each datacase available in the volatile catalogue
        for datacase in self._getsortedcataloguekeys(cataloguekeyslist=list(volcat.keys())):
            mechanisms = volcat[datacase].keys()

            # Extend the dictionary with a key for each data case
            znd_pratio_dict[datacase] = {}

            # For each mechanism in a datacase
            for mechanism in mechanisms:

                # Work out the quantities required by the Shock and Detonation Toolbox
                p1, t1, q, diluentmoles = self.case_initialconditions(datacase=datacase, mechanism=mechanism)

                # If 'fast' is true and q=C2H4:3O2:xAr, look for a pre-baked approximation
                if (fast is True) and ("C2H4:1.00 O2:3.00 Ar:" in q) and ("Ar" == q.split(" ")[-1].split(":")[0]):

                    # T = 288 K, P = avg(10-100 kPa) ~ 30 kPa
                    a1 = -0.130574607
                    b1 = 1.770783750

                    # T = 298 K, P = 100 k Pa
                    a2 = -0.122975381
                    b2 = 1.783500000

                    # Compute weighting, to decide which method to use
                    p_weight = np.interp(p1, [30000, 100000], [0, 1])
                    t_weight = np.interp(t1, [288, 298], [0, 1])
                    pt_weight = p_weight + t_weight / 2

                    a = a1 * (1 - pt_weight) + a2 * pt_weight
                    b = b1 * (1 - pt_weight) + b2 * pt_weight
                    vncjratio = a * (math.log(diluentmoles, 10) - 1) + b

                    # Calculate values from scratch (very slow!)
                else:
                    # Set up an object of the Shock & Detonation Toolbox
                    zerod_obj = SDT0d(p_initial=p1, t_initial=t1, q_initial=q)

                    # Find the pressure ratios (this takes a long time because CJ speed takes ages to calculate)
                    vncjratio = zerod_obj.calc_vnpressure() / zerod_obj.calc_cjpressure()

                # Store the vncj ratio into the dictionary
                znd_pratio_dict[datacase][mechanism] = vncjratio

        # If a SDT0d object was created, destroy it
        try:
            del zerod_obj
        except UnboundLocalError:
            pass

        return znd_pratio_dict


if __name__ == "__main__":
    test1d = AMROC1d()

    moles = [10]
    for mole in moles:
        test1d.case_read(datacasefolder=f"50_10_{mole}_data")
        # test1d.case_read(datacasefolder=f"90_50_{mole}_data")

    test1d.data_detonationcheck()

    # for testcase in test1d.datacatalogue:
    #     test1d.case_read(datacasefolder=testcase)
    #
    # test1d.export_statistics()
    # test1d.export_detonationreport()
    # test1d.plot_mechanisms(colourseed=True)
