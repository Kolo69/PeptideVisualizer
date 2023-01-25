#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (c) 2022
#     Matej Kolarič (matej.kolaric@ijs.si), Robert Vidmar, Marko Fonović
#     Biochemistry, molecular and structural biology - B1, Jozef Stefan Institute, Jamova cesta 39, Ljubljana Slovenia
#     under GPLv3 (http://www.gnu.org/licenses/gpl-3.0.html)
#
# BETA testers: Tilen Sever, Marija Grozdanić, Sara Ivanovski, Andreja Kozak
#
#

##############################################################################
#               Nastavitve
##############################################################################

DEBUG = False
NOTICE = False
SELFUPDATE = True
GUI = True
DARK_MODE = False
bar_height = 0.9
bar_width = 0.8
alpha_boost = 0.05
default_peptide_threshold = 5 # v odstotkih
uniProtFeature_coeficient = 0.07 # +x% višine grafa je legenda
UniProt_Features_alpha = 0.5
replicates_MAX = 10000
transparent_bg = False
TXT_OUTPUT = True
Multiproc = True
N_download_threads = 16 # včasih zlaufa lepo z 100+, ampak obstaja verjetnost da se sesuje
N_write_threads = 1 # Nima bistvenega vpliva + pri malem številu zaznanih proteinov včasih neki zjebe
DEV_DPI = 72.0
DEV_scale = 1.0

UniProt_Feature_types = ["Processing", "Region", "Secnondary structure"]
UniProt_Processing_Features = ["Signal peptide", "Transit peptide", "Propeptide", "Chain", "Peptide"]
UniProt_Region_Features = ["Topological domain", "Transmembrane", "Intramembrane", "Domain", "Repeat", "Calcium binding", "Zinc finger", "DNA binding", "Nucleotide binding", "Region", "Coiled coil", "Motif", "Compositional bias"]
UniProt_2ndary_Features = ["Helix", "Turn", "Beta strand"]

LEGEND_Processing = [
    ['Signal peptide', '#ff0000'],
    ['Transit peptide', '#FFA500'],
    ['Propeptide', '#FFA500'],
    ['Chain', '#00ff00'],
    ['Peptide', '#FFA500']
]

LEGEND_Region = [
    ['Topological domain', '#00aa00'],
    ['Domain', '#00aa00'],
    ['Region', '#ffff00'],
    ['Repeat', '#00ffff'],
    ['Calcium binding', '#aa0000'], 
    ['Zinc finger', '#aa0000'], 
    ['DNA binding', '#aa0000'],
    ['Nucleotide binding', '#ff0000'],
    ['Intramembrane', '#0000ff'],
    ['Transmembrane', '#0000ff'],
    ['Coiled coil', '#0000ff'],
    ['Motif', '#00aa00'],
    ['Compositional bias', '#ff0000']
]

LEGEND_2ndary = [
    ['Helix', '#ff0000'],
    ['Turn', '#00ff00'],
    ['Beta strand', '#0000ff']
]

LEGEND = [
    #hidden
    ['Uncolored', '#ffffff'],
    ['Alternative sequence', '#ffffff']
] + LEGEND_Processing + LEGEND_Region + LEGEND_2ndary

MW_PRESET = {
    '12 % Tris-Glycine': {
        5: [100, 55, 35, 20],
        6: [120, 60, 40, 25, 15],
        7: [130, 70, 45, 35, 25, 15],
        8: [140, 80, 50, 40, 30, 20, 15],
        9: [150, 85, 60, 40, 35, 25, 20, 10],
        10: [170, 95, 70, 50, 40, 30, 25, 15, 10],
        11: [185, 100, 70, 50, 40, 35, 35, 20, 15, 10],
        12: [200, 110, 75, 60, 45, 35, 32, 25, 20, 15, 10],
        13: [210, 120, 85, 60, 50, 40, 35, 30, 25, 20, 15, 10],
        14: [220, 130, 90, 70, 55, 45, 40, 35, 25, 22, 20, 15, 10],
        15: [230, 135, 100, 75, 60, 50, 40, 35, 30, 25, 22, 18, 15, 8],
        16: [240, 140, 100, 75, 60, 50, 45, 40, 35, 28, 25, 21, 15, 12, 7]
    },
    '8-16 % Tris-Glycine': {
        5: [100, 45, 25, 10],
        6: [120, 60, 35, 20, 10],
        7: [140, 70, 45, 28, 18, 8],
        8: [160, 80, 50, 35, 23, 15, 7],
        9: [175, 90, 60, 40, 30, 20, 13, 7],
        10: [190, 100, 65, 45, 35, 25, 18, 12, 6],
        11: [200, 110, 75, 55, 40, 30, 23, 15, 10, 6],
        12: [210, 120, 80, 60, 45, 35, 27, 20, 15, 9, 5],
        13: [225, 130, 90, 65, 50, 40, 32, 25, 20, 13, 8, 5],
        14: [230, 140, 100, 70, 55, 45, 35, 28, 22, 17, 12, 8, 4],
        15: [240, 150, 105, 75, 60, 47, 40, 32, 25, 20, 15, 11, 8, 4],
        16: [250, 160, 110, 80, 60, 50, 42, 35, 28, 23, 20, 15, 10, 7, 4]
    }
}

PROGRAM_NAME = 'PeptideVisualizer'
VERSION = '1.2c'

online_srv = ["8.8.8.8", "9.9.9.9", "208.67.222.222", "1.1.1.1"] # Google, Quad9, OpenDNS Home, Cloudflare
update_srv = "proteom.eu"
update_url = f"https://proteom.eu/{PROGRAM_NAME}/{PROGRAM_NAME}.py"

README = f"""##########################################################################################################

{PROGRAM_NAME} v{VERSION}

(c) 2022, Matej Kolarič (matej.kolaric@ijs.si), Robert Vidmar, Marko Fonović
Jozef Stefan Institute, Slovenia

PeptideVisualizer is a tool for visualization of proteins in proteomic datasets based on the coverage of mass spectrometry identified peptides. The PeptideVisualizer principle originates from the PROTOMAP (Dix MM, Simon GM, Cravatt BF. Global mapping of the topography and magnitude of proteolytic events in apoptosis. Cell. 2008) algorithm with an addition of several upgrades as:
  - enables visualization of three data dimensions (x: protein sequence coverage, y: position in each fraction, z: intensity based colour density)
  - enables input data from peptides.txt table generated by MaxQuant (Tyanova S, Temu T, Cox J. The MaxQuant computational platform for mass spectrometry-based shotgun proteomics. Nat Protoc. 2016)
  - enables either Intensity or LFQ intensity values input
  - generates Uniprot data mined protein domains and motifs for improved visualization of curated Uniprot protein features. 

##########################################################################################################"""

#globalne spremenljivke
usr_name = ""
DataGroup1_name = ""
DataGroup2_name = ""
N_max = 0
replicates = 1
error_bar = True
horizontal_lines = True
SVG_OUTPUT = True
use_LFQ = False
useUniProt = False
setMW = False
group_normalization = False
useUniProt_features = []
DataGroup1 = []
DataGroup2 = []
data_indexes = []
data_titles = []
setMW_Labels = []
peptide_threshold = default_peptide_threshold
DEBUG_txt = []

empty_char = ' '
start_char = '├'
mid_char = '-'
end_char = '┤'
cross_char = '┼'

##############################################################################
#               Moduli
##############################################################################
print("\n"+README+"\n")

DEBUG_txt.append(f"Initilizing {PROGRAM_NAME} v{VERSION} ...")
if NOTICE:
    print(DEBUG_txt[-1])

# basic stuff
import os
os.system("") #colorfix
from os.path import exists
import shutil
import csv
import subprocess
import sys
from sys import platform
import re
import time
import calendar
from datetime import datetime
import unicodedata
import codecs
import json
import math
import statistics
from multiprocessing.pool import ThreadPool

# output formating + pause in case of missing lib
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

if platform == "linux" or platform == "linux2":
    OS = "Linux"
elif platform == "darwin":
    OS = "MacOS"
elif platform == "win32":
    OS = "Windows"
else:
    OS = "Other"
DEBUG_txt.append(f"NOTICE: OS = {OS}")
if NOTICE:
    print(DEBUG_txt[-1])

def pause():
    DEBUG_txt.append(f"Waiting for user to acknowledge ...")
    programPause = input("Press the <ENTER> key to continue...")

### DEPENDENCIES REQUIRED FOR UPGRADE ####

# Preverit če zlaufa "pip -V", samo v tem primeru bo namestilo
DEBUG_txt.append('NOTICE: Checking if pip is installed ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    subprocess.check_call([sys.executable, "-m", "pip", "-V", "--quiet"],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
except:
    DEBUG_txt.append("ERROR! Python module pip not found!")
    print(f"{bcolors.FAIL}{DEBUG_txt[-1]} {PROGRAM_NAME} is unable to install it You will have to install it manually.{bcolors.ENDC}")
    if OS == "Windows":
        print(f"Please run \"sudo apt-get install python3-pip\" or equivalent to install python3-pip package!\n")
    else:
        print(f"I have no idea how pip is not installed on OS={OS} ... Please try updating python.")
    pause()
    exit()

# REQUIRED: numpy
DEBUG_txt.append('NOTICE: Importing numpy ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    import numpy as np
except ImportError:
    DEBUG_txt.append('WARNING! Python module numpy not found!')
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install numpy\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling numpy ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'numpy', "--quiet"])
    time.sleep(1)
    #import numpy as np # zašteka včasih!
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed numpy! Please restart the script.")
    print(DEBUG_txt[-1]+"\n")
    pause()
    exit()

# REQUIRED: matplotlib
DEBUG_txt.append('NOTICE: Importing matplotlib ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    import matplotlib
except ImportError:
    DEBUG_txt.append('WARNING! Python module matplotlib not found!')
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install matplotlib\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling matplotlib ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'matplotlib', "--quiet"])
    time.sleep(1)
    #import matplotlib
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed matplotlib! Please restart the script")
    print(DEBUG_txt[-1]+"\n")
    pause()
    exit()

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
matplotlib.use("Agg") #brez te vrstice mi odleti vn iz memory limita (2GB). Agg, is a non-interactive backend that can only write to files.

# REQUIRED: requests - internetih čoj!
DEBUG_txt.append('NOTICE: Importing requests ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    import requests
except ImportError:
    DEBUG_txt.append('WARNING! Python module requests not found!')
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install requests\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling requests ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'requests', "--quiet"])
    time.sleep(1)
    #import requests
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed requests! Please restart the script")
    print(DEBUG_txt[-1]+"\n")
    pause()
    exit()

# datumi za apdejt
DEBUG_txt.append('NOTICE: Importing dateutil.parser ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    from dateutil.parser import parse as parsedate
except ImportError:
    DEBUG_txt.append("WARNING! Python module dateutil not found!")
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install python-dateutil\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling python-dateutil ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'python-dateutil', "--quiet"])
    time.sleep(1)
    #from dateutil.parser import parse as parsedate
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed dateutil! Please restart the script.")
    print(DEBUG_txt[-1]+"\n")
    pause()
    exit()

# REQUIRED: alive_progress
DEBUG_txt.append('NOTICE: Importing alive_progress ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    from alive_progress import alive_bar
except ImportError:
    DEBUG_txt.append("WARNING! Python module alive-progress not found!")
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install alive-progress\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling alive-progress ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'alive-progress', "--quiet"])
    time.sleep(1)
    #from alive_progress import alive_bar
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed alive_progress! Please restart the script.")
    print(DEBUG_txt[-1]+"\n")
    pause()
    exit()

###
# GUI (tkinter) import + basic checks
###
if GUI:
    DEBUG_txt.append('NOTICE: Importing FigureCanvasTkAgg ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #for MW setting
    except:
        DEBUG_txt.append("WARNING! Python module python3-pil.imagetk not found!")
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install pillow\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
        pause()
        print("Istalling alive-progress ... ")
        subprocess.check_call([sys.executable, "-m", "pip", "install", 'pillow', "--quiet"])
        time.sleep(1)
        #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed python3-pil.imagetk! Please restart the script.")
        print(DEBUG_txt[-1]+"\n")
        pause()
        exit()
    
    DEBUG_txt.append('NOTICE: Importing tkinter ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        import tkinter
    except ImportError:
        GUI = False
        DEBUG_txt.append("WARNING! Python module tkinter (python3-tk) not found!")
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} is unable to install it for you.{bcolors.ENDC}")
        if OS == "Linux":
            print(f"Please run \"sudo apt-get install python3-tk\" or equivalent to install python3-tk ...")
        elif OS == "MacOS":
            print(f"Please run \"brew install python-tk\" to install python3-tk ...")
        else:
            print(f"I have no idea how tkinter is not installed on OS={OS} ... Please try updating python.")
        pause()
        exit(2)
    
    DEBUG_txt.append('NOTICE: Importing askopenfilename ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        from tkinter.filedialog import askopenfilename
        #root = tkinter.Tk()
        #root.withdraw() # hack to hide tkinter window NOTE: puts windows out of focus
        #root.after(10, lambda: root.withdraw())
    except:
        GUI = False
if GUI:
    try:
        import tkinter as tk
    except ImportError:
        import Tkinter as tk

    import ctypes #lahko to kje sfaila?
    if OS == "Windows":
        try:
            ctypes.windll.shcore.SetProcessDpiAwareness(2) # win version >= 8.1
            DEBUG_txt.append(f'NOTICE: windll.shcore.SetProcessDpiAwareness(2); Windows ver. >= 8.1')
        except:
            ctypes.windll.user32.SetProcessDPIAware() # win 8.0 or less
            DEBUG_txt.append(f'NOTICE: windll.user32.SetProcessDPIAware(); Windows ver. <= 8.0') #kako se to obnaša na drugih OSjih??
        if NOTICE:
            print(DEBUG_txt[-1])

    try:
        tk.Tk().destroy()
    except:
        GUI = False

GUI_issue = False
if GUI:
    tk_ver = tk.Tcl().call("info", "patchlevel").split(".")
    DEBUG_txt.append(f'NOTICE: tkinter version {".".join(tk_ver)} deteccted.')
    if NOTICE:
        print(DEBUG_txt[-1])
    if int(tk_ver[0]) == 8 and int(tk_ver[1]) == 6 and int(tk_ver[2]) <= 11:
        GUI_issue = True
        DEBUG_txt.append(f'{bcolors.WARNING}WARNING: {PROGRAM_NAME} might expirience display difficulties with tkinter v.{".".join(tk_ver)}. Try updating python/tkinter.{bcolors.ENDC}')
        print(DEBUG_txt[-1])

DEBUG_txt.append(f'NOTICE: Graphical User Interface = {GUI}')
if NOTICE:
    print(DEBUG_txt[-1])

##############################################################################
#               Funkcije
##############################################################################


def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value)
    return re.sub(r'[-\s]+', '-', value).strip('-_')

def sequence_adjust(x, pos):
    if pos == 'mid':
        if x == cross_char :
            return cross_char
        elif x == empty_char or x == mid_char :
            return mid_char
        else:
            return x
    elif pos == 'start':
        if x == empty_char or x == mid_char or x == start_char :
            return start_char
        else:
            return cross_char
    elif pos == 'end':
        if x == empty_char or x == mid_char or x == end_char :
            return end_char
        else:
            return cross_char
    else:
        return 'E'

def replace_range_by_ints(m):
    a = m.group(1)
    b = m.group(2)
    return ','.join(str(i) for i in range(int(a), int(b) + 1))

def userUniProtFeatures():
    user_features = []
    pattern = re.compile('(\d+)\-(\d+)')
    print(f"{bcolors.OKBLUE}\t*** Available features: ***\nID\t=>NAME{bcolors.ENDC}")
    for feature_id in range(len(UniProt_Feature_types)):
        print (feature_id+1, "\t=> ", UniProt_Feature_types[feature_id])
    
    DEBUG_txt.append(f"NOTICE: Waiting for user input to form integer list (min: 0, max: {len(UniProt_Feature_types)})")
    if NOTICE: 
        print(f"{DEBUG_txt[-1]}\n")
    else:
        print("")
    usr_in = input("Please enter comma separated IDs of features you want to include [default: none]: ").replace(" ", "").replace("\t", "").replace("\n", "")
    if not usr_in:
        return False
    else:
        val = re.sub(pattern, replace_range_by_ints, usr_in).split(",") # ločim vejice in spremenim - v seznam
        for single in val:
            try:
                int(single)
            except ValueError:
                print(f"{bcolors.FAIL}ERROR! Not a valid input! Terminating ...{bcolors.ENDC}")
                pause()
                exit(2)
            if int(single)-1 < len(UniProt_Feature_types):
                user_features.append(int(single)-1)
    user_features.sort()
    DEBUG_txt.append("NOTICE: Using features ")
    for feature_id in user_features:
        DEBUG_txt[-1] += f"\"{UniProt_Feature_types[feature_id]}\", "
    DEBUG_txt[-1] = DEBUG_txt[-1][:-2]
    print(f"{bcolors.OKBLUE}\t{DEBUG_txt[-1]}{bcolors.ENDC}")
    
    return(user_features)

def userGrouping(data):    
    global DataGroup1_name
    global DataGroup2_name
    print(f"{bcolors.OKBLUE}\t*** Detected experiments: ***\nID\t=>NAME\n{bcolors.ENDC}")
    for exp_key in range(len(data)):
        print (exp_key+1, "\t=> ", data[exp_key])

    pattern = re.compile('(\d+)\-(\d+)') # za spreminjanje 1-5 v 1,2,3,4,5
    
    print(f"\n{bcolors.OKBLUE}INPUT EXAMPLES:\t1,2,3\t\t1-3\t\t1-5,7,8-10,11\t...{bcolors.ENDC}\n")
    
    ##### DataGroup1 #####
    usr_in = input("DataGroup1 comma separated IDs list: ").replace(" ", "").replace("\t", "").replace("\n", "")
    if not usr_in:
        print(f"{bcolors.FAIL}ERROR! No user input! Terminating ...{bcolors.ENDC}")
        pause()
        exit(2)

    val = re.sub(pattern, replace_range_by_ints, usr_in).split(",") # ločim vejice in spremenim - v seznam
    for single in val:
        try:
            data_id = int(single)-1
        except ValueError:
            print(f"{bcolors.FAIL}ERROR! Not a valid input! Terminating ...{bcolors.ENDC}")
            pause()
            exit(2)
        
        if data_id >= 0 and data_id <= len(data):
            if data_id not in DataGroup1:
                DataGroup1.append(data_id)
            else:
                DEBUG_txt.append(f"WARNING: Data entry \"'{data_id+1}'\" is entered twice! Removing second instance!")
                print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
        else:
            DEBUG_txt.append(f"WARNING: Data entry \"{data_id+1}\" is out of range! Removing it!")
            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")

    DEBUG_txt.append(f"Entry successful! DataGroup1: ")
    for key in DataGroup1:
        DEBUG_txt[-1] += data[key]+", "
    DEBUG_txt[-1] = DEBUG_txt[-1][:-2]
    print(f"{bcolors.OKBLUE}{DEBUG_txt[-1]}{bcolors.ENDC}")
    
    DataGroup1_name = slugify(input("Enter DataGroup1 name: "))
    DEBUG_txt.append(f"NOTICE: DataGroup1_name='{DataGroup1_name}'")
    if NOTICE:
        print(DEBUG_txt[-1])
    
    ##### DataGroup2 #####
    print("")
    usr_in = input(f"DataGroup2 comma separated IDs list (0 or {len(DataGroup1)} elements required): ").replace(" ", "").replace("\t", "").replace("\n", "")
    if usr_in:
        val = re.sub(pattern, replace_range_by_ints, usr_in).split(",")
        for single in val:
            try:
                data_id = int(single) - 1
            except ValueError:
                print(f"{bcolors.FAIL}ERROR! Not a valid input! Terminating ...{bcolors.ENDC}")
                pause()
                exit(2)
            
            if data_id >= 0 and data_id < len(data):
                if data_id not in DataGroup2:
                    DataGroup2.append(data_id)
                else:
                    DEBUG_txt.append(f"WARNING: Data entry \"{data_id+1}\" is entered twice! Removing second instance!")
                    print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
            else:
                DEBUG_txt.append(f"WARNING: Data entry \"{data_id+1}\" is out of range! Removing it!")
                print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")

        #print(DataGroup2)
        if len(DataGroup1) != len(DataGroup2):
            print(f"{bcolors.FAIL}ERROR! DataGroups not the same size! Terminating ...{bcolors.ENDC}")
            pause()
            exit(2)
        
        DEBUG_txt.append(f"Entry successful! DataGroup2: ")
        for key in DataGroup2:
            DEBUG_txt[-1] += data[key]+", "
        DEBUG_txt[-1] = DEBUG_txt[-1][:-2]
        print(f"{bcolors.OKBLUE}{DEBUG_txt[-1]}{bcolors.ENDC}")
        
        DataGroup2_name = slugify(input("Enter DataGroup2 name: "))
        DEBUG_txt.append(f"NOTICE! DataGroup2_name='{DataGroup2_name}'")
        if NOTICE:
            print(DEBUG_txt[-1])

    return DataGroup1, DataGroup2

def query_yes_no(question, default="yes"):
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

###
# GUI v1.2
###
class PeptideVisualizerGUI(tk.Tk):
    def __init__(self, *args, **kwargs):
        global ONLINE
        tk.Tk.__init__(self, *args, **kwargs)
        
        # DARK/LIGHT Mode
        if DARK_MODE:
            self.bg, self.fg, self.ac1, self.ac2 = ('#282828', 'white', '#404040', '#B3B3B3')
        else:
            self.bg, self.fg, self.ac1, self.ac2 = ('#f0f0ed', 'black', '#ffffff', '#E9DAC1')
        
        #self.config(bg="#26242f")
        
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()
        screen_dpi = self.winfo_fpixels('1i')
        OS_scale = screen_dpi / 96 # 96 is a scale factor for 100%
        scale_f = screen_dpi / DEV_DPI
        screen_w_eff = round(screen_w / OS_scale)
        screen_h_eff = round(screen_h / OS_scale)
        wh = round(650 * OS_scale)
        ww = round(600 * OS_scale)
        
        DEBUG_txt.append(f"NOTICE: resolution={screen_w}x{screen_h} dpi={screen_dpi} OS_scale={OS_scale} ({screen_w_eff}x{screen_h_eff})")
        DEBUG_txt.append(f"NOTICE: scale_f={scale_f} => {wh}x{ww}")
        if NOTICE:
            print(DEBUG_txt[-2]+"\n"+DEBUG_txt[-1])
        
        self.tk.call('tk', 'scaling', scale_f)
        self.font = ('Serif', round(12*scale_f))
        self.configure(bg=self.bg, borderwidth=0)
        self.title('PeptideVisualizer')
        #self.geometry(f"{wh}x{ww}")
        self.resizable()
        #self.resizable(width=False, height=False)
        
        # Variables
        self.Experiment_name = tk.StringVar()
        self.replicates_var = tk.StringVar()
        self.replicates_var.set(1)
        self.peptide_threshold_var = tk.StringVar()
        self.peptide_threshold_var.set(default_peptide_threshold)
        self.saveSVG = tk.BooleanVar(value=True)
        
        self.DataSets_vars = tk.StringVar()
        
        self.DataGroup1_name = tk.StringVar()
        self.DataGroup1_name.set("DataGroup 1")
        self.DataGroup1_vars = tk.StringVar()
        
        self.DataGroup2_name = tk.StringVar()
        self.DataGroup2_name.set("DataGroup 2")
        self.DataGroup2_vars = tk.StringVar()
        
        self.useLFQ = tk.BooleanVar()
        
        self.error_bar = tk.BooleanVar(value=True)
        self.horizontal_lines = tk.BooleanVar(value=True)
        self.setMW = tk.BooleanVar(value=True)
        self.group_normalization = tk.BooleanVar()
        
        if ONLINE:
            self.useUniProt = tk.BooleanVar(value=True)
        else:
            self.useUniProt = tk.BooleanVar()
        
        self.useFeatures1 = tk.BooleanVar(value=True)
        self.useFeatures2 = tk.BooleanVar()
        self.useFeatures3 = tk.BooleanVar(value=True)
        
        # Layout - Experiment_label
        self.Experiment_label = tk.Label(text="Experiment Name:", width=32, anchor="e")
        self.Experiment_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.Experiment_label.grid(row=0, column=0)
        
        self.Experiment_text = tk.Entry(textvariable = self.Experiment_name,width=32)
        self.Experiment_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.Experiment_text.grid(row=0, column=1, pady=(10, 10))
        
        self.DataSets_name = tk.Label(text="Data from peptides.txt", width=32, anchor="w")
        self.DataSets_name.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.DataSets_name.grid(row=1, column=0)
        
        self.DataGroup1_text = tk.Entry(textvariable = self.DataGroup1_name,width=32)
        self.DataGroup1_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup1_text.grid(row=1, column=1)
        
        self.DataGroup2_text = tk.Entry(textvariable = self.DataGroup2_name,width=32)
        self.DataGroup2_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup2_text.grid(row=1, column=2)
        
        # DataSets retrived from peptides.txt
        self.DataSets_list = tk.Listbox(listvariable=self.DataSets_vars, height=21, width=32, selectmode='extended')
        self.DataSets_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataSets_list.grid(row=2, column=0)
        self.DataSets_vars.set(value=data_titles_INT)
        
        # DataGroup1
        self.DataGroup1_list = tk.Listbox(listvariable=self.DataGroup1_vars, height=21, width=32, selectmode='extended')
        self.DataGroup1_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup1_list.grid(row=2, column=1)
        self.add1 = tk.Button(text=' ADD ', command=self.add2one)
        self.add1.configure(bg=self.ac1, fg=self.fg)
        self.add1.grid(row=3, column=1, padx=(0,50))
        self.remove1 = tk.Button(text=' DEL ', command=self.del_one)
        self.remove1.configure(bg=self.ac1, fg=self.fg)
        self.remove1.grid(row=3, column=1, padx=(50,0))
        
        # DataGroup2
        self.DataGroup2_list = tk.Listbox(listvariable=self.DataGroup2_vars, height=21, width=32, selectmode='extended')
        self.DataGroup2_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup2_list.grid(row=2, column=2)
        self.add2 = tk.Button(text=' ADD ', command=self.add2two)
        self.add2.configure(bg=self.ac1, fg=self.fg)
        self.add2.grid(row=3, column=2, padx=(0,50))
        self.remove2 = tk.Button(text=' DEL ', command=self.del_two)
        self.remove2.configure(bg=self.ac1, fg=self.fg)
        self.remove2.grid(row=3, column=2, padx=(50,0))
        
        # Everything else - first column
        self.useLFQ_check = tk.Checkbutton(text="use LFQ", variable=self.useLFQ, width=28, anchor="w", onvalue=True, offvalue=False)
        self.useLFQ_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useLFQ_check.grid(row=4, column=0, pady=(10,0))
        
        self.replicates_label = tk.Label(text="Replicates   N = ", anchor="w", width=30)
        self.replicates_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.replicates_label.grid(row=5, column=0)
        self.replicates_var_text = tk.Entry(textvariable = self.replicates_var, width=3)
        self.replicates_var_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.replicates_var_text.grid(row=5, column=0)
        
        self.error_bar_check = tk.Checkbutton(text="Display error bars", variable=self.error_bar, width=28, anchor="w", onvalue=True, offvalue=False, state="disabled")
        self.error_bar_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.error_bar_check.grid(row=6, column=0)
        
        self.peptide_threshold_label = tk.Label(text="Threshold    T =        %", anchor="w", width=30)
        self.peptide_threshold_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.peptide_threshold_label.grid(row=6, column=0)
        self.peptide_threshold_var_text = tk.Entry(textvariable = self.peptide_threshold_var, width=3)
        self.peptide_threshold_var_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.peptide_threshold_var_text.grid(row=6, column=0)
        
        self.saveSVG_check = tk.Checkbutton(text="Save as .svg", variable=self.saveSVG, width=28, anchor="w", onvalue=True, offvalue=False)
        self.saveSVG_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.saveSVG_check.grid(row=7, column=0)
        
        # second column
        self.horizontal_lines_check = tk.Checkbutton(text="Lines between slices", variable=self.horizontal_lines, width=28, anchor="w", onvalue=True, offvalue=False)
        self.horizontal_lines_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.horizontal_lines_check.grid(row=4, column=1, pady=(10,0))
        
        self.setMW_check = tk.Checkbutton(text="Molar mass on Y-ax", variable=self.setMW, width=28, anchor="w", onvalue=True, offvalue=False)
        self.setMW_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.setMW_check.grid(row=5, column=1)
        
        self.group_normalization_check = tk.Checkbutton(text="Group normalization", variable=self.group_normalization, width=28, anchor="w", onvalue=True, offvalue=False)
        self.group_normalization_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.group_normalization_check.grid(row=6, column=1)
        
        #third column
        self.useUniProt_check = tk.Checkbutton(text="use UniProt online DB", variable=self.useUniProt, width=28, anchor="w", onvalue=True, offvalue=False, command=self.chgUniProt_check)
        self.useUniProt_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures1_check = tk.Checkbutton(text="show "+UniProt_Feature_types[0], variable=self.useFeatures1, width=28, anchor='w', onvalue=True, offvalue=False)
        self.useFeatures1_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures2_check = tk.Checkbutton(text="show "+UniProt_Feature_types[1], variable=self.useFeatures2, width=28, anchor="w", onvalue=True, offvalue=False)
        self.useFeatures2_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures3_check = tk.Checkbutton(text="show "+UniProt_Feature_types[2], variable=self.useFeatures3, width=28, anchor="w", onvalue=True, offvalue=False)
        self.useFeatures3_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        
        if not ONLINE:
            self.useUniProt_check.config(text="UniProt unavailable (you are offline)", state="disabled")
            self.useFeatures1_check.config(state="disabled")
            self.useFeatures2_check.config(state="disabled")
            self.useFeatures3_check.config(state="disabled")  
        self.useUniProt_check.grid(row=4, column=2, pady=(10,0))
        self.useFeatures1_check.grid(row=5, column=2, padx=(10,0))
        self.useFeatures2_check.grid(row=6, column=2, padx=(10,0))
        self.useFeatures3_check.grid(row=7, column=2, padx=(10,0))
        
        self.VisualizerDo_btn = tk.Button(text="Visualize", command=self.VisualizerDo)
        self.VisualizerDo_btn.configure(bg=self.ac1, fg=self.fg)
        self.VisualizerDo_btn.grid(row=15, column=1, pady=(10,0))
    
    def chgUniProt_check(self):
        if self.useUniProt.get() == True:
            self.useFeatures1_check.config(state="normal")
            self.useFeatures2_check.config(state="normal")
            self.useFeatures3_check.config(state="normal")
        else:
            self.useFeatures1_check.config(state="disabled")
            self.useFeatures2_check.config(state="disabled")
            self.useFeatures3_check.config(state="disabled")
    
    def add2one(self):
        if self.DataSets_list.curselection() == ():
            return

        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup1_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        #result_list = sorted(list(set(existing_list + selected_list)))
        
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup1_vars.set(value=result_list)
        
    def del_one(self):
        if self.DataGroup1_list.curselection() == ():
            return
        selection = self.DataGroup1_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup1_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup1_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup1_vars.set(value=existing_list)
    
    def add2two(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup2_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup2_vars.set(value=result_list)
        
    def del_two(self):
        if self.DataGroup2_list.curselection() == ():
            return
        selection = self.DataGroup2_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup2_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup2_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup2_vars.set(value=existing_list)
    
    def VisualizerDo(self):
        global usr_name
        global DataGroup1_name
        global DataGroup2_name
        global N_max
        global replicates
        global error_bar
        global horizontal_lines
        global SVG_OUTPUT
        global use_LFQ
        global useUniProt
        global setMW
        global group_normalization
        global useUniProt_features
        global DataGroup1
        global DataGroup2
        global data_indexes
        global data_titles
        global setMW_Labels
        global peptide_threshold
        
        DEBUG_txt.append(f'NOTICE: Acquiring user settings from GUI...')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        # Name
        usr_name = slugify(self.Experiment_name.get())
        
        DEBUG_txt.append(f'NOTICE: usr_name = "{usr_name}"')
        if NOTICE:
            print(DEBUG_txt[-1])

        #LFQ
        use_LFQ = self.useLFQ.get()
        if existsLFQ and use_LFQ:
            data_indexes = data_indexes_LFQ
            data_titles = data_titles_LFQ
        else:
            data_indexes = data_indexes_INT
            data_titles = data_titles_INT
            
        DEBUG_txt.append(f'NOTICE: use_LFQ = {use_LFQ}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        #DataGroups
        DataGroup1_name = self.DataGroup1_name.get()
        DataGroup2_name = self.DataGroup2_name.get()
        
        DataGroup1_names = [line.strip(' \'') for line in self.DataGroup1_vars.get()[1:-1].split(',')]
        while ("" in DataGroup1_names):
            DataGroup1_names.remove("")
        DataGroup2_names = [line.strip(' \'') for line in self.DataGroup2_vars.get()[1:-1].split(',')]
        while ("" in DataGroup2_names):
            DataGroup2_names.remove("")
            
        for t in DataGroup1_names:
            #DataGroup1.append(data_indexes[data_titles.index(t)])
            DataGroup1.append(data_titles.index(t))
        for t in DataGroup2_names:
            DataGroup2.append(data_titles.index(t))
        
        DEBUG_txt.append(f"NOTICE: {DataGroup1_name} = {json.dumps(DataGroup1_names)}\nNOTICE: {DataGroup2_name} = {json.dumps(DataGroup2_names)}") #v prejšnjem debugu izpisujem DataGroup2 in ne DataGroup2_names
        if NOTICE:
            print(DEBUG_txt[-1])
        
        if len(DataGroup1) >= 1:
            N_max = len(DataGroup1)
            if N_max != len(DataGroup2):
                DataGroup2 = []
                print(f"{bcolors.WARNING}WARNING! Only DataGroup1 is defined!{bcolors.ENDC}")
        else:
            print(f"{bcolors.FAIL}ERROR! Not a single group is defined!{bcolors.ENDC}")
            return
        # Replicates
        try:
            replicates = int(self.replicates_var.get())
        except ValueError:
            print(f"{bcolors.WARNING}WARNING! Number of replicates could not be interpreted. Using default value 1{bcolors.ENDC}")
            replicates = 1
        if replicates < 1 or replicates > replicates_MAX:
            print(f"{bcolors.WARNING}WARNING! Number of replicates makes no sense ({replicates}). Using default value 1{bcolors.ENDC}")
            replicates = 1
        if N_max % replicates != 0:
            print(f"{bcolors.FAIL}ERROR! Number of replicates ({replicates}) makes no sense considering group size ({N_max}). {N_max} is not divisible by {replicates}!{bcolors.ENDC}")
            return
        DEBUG_txt.append(f'NOTICE: replicates = {replicates}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        N_max = int(N_max / replicates) # nov maximum
            
        #Threshold
        try:
            peptide_threshold = float(self.peptide_threshold_var.get())
        except ValueError:
            DEBUG_txt.append(f"WARNING! Relative peptide detection threshold could not be interpreted. Using default value {default_peptide_threshold}%")
            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
            peptide_threshold = default_peptide_threshold
        if peptide_threshold < 0 or peptide_threshold > 100:
            DEBUG_txt.append(f"WARNING! Relative peptide detection threshold makes no sense ({peptide_threshold}%). Using default value {default_peptide_threshold}%")
            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
            peptide_threshold = default_peptide_threshold

        DEBUG_txt.append(f'NOTICE: peptide_threshold = {peptide_threshold}%')
        if NOTICE:
            print(DEBUG_txt[-1])

        useUniProt = self.useUniProt.get()
        if useUniProt:
            if self.useFeatures1.get():
                useUniProt_features.append(0)
            if self.useFeatures2.get():
                useUniProt_features.append(1)
            if self.useFeatures3.get():
                useUniProt_features.append(2)
        
        DEBUG_txt.append(f'NOTICE: useUniProt = {useUniProt}; useUniProt_features = {json.dumps(useUniProt_features)}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        SVG_OUTPUT = self.saveSVG.get()
        DEBUG_txt.append(f'NOTICE: SVG_OUTPUT = {SVG_OUTPUT}')
        if NOTICE:
            print(DEBUG_txt[-1])
            
        error_bar = self.error_bar.get()
        DEBUG_txt.append(f'NOTICE: error_bar = {error_bar}')
        if NOTICE:
            print(DEBUG_txt[-1])
            
        horizontal_lines = self.horizontal_lines.get()
        DEBUG_txt.append(f'NOTICE: horizontal_lines = {horizontal_lines}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        self.quit()
        if OS == "Windows":
            self.destroy() # uničim okno, da ne bosta dve okni eno čez drugo in če ne uničiš ma windows fatal error
        
        setMW = self.setMW.get()
        if setMW and N_max > 1:
            App_setMW = PeptideVisualizerSetMw()
            App_setMW.mainloop()
            if setMW_Labels == []:
                setMW = False
                DEBUG_txt.append(f"WARNING! No Y-ax labels given. Y-ax will be labeled by slices. setMW = {setMW}")
                print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
            
        DEBUG_txt.append(f'NOTICE: setMW = {setMW}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        group_normalization = self.group_normalization.get()
        DEBUG_txt.append(f'NOTICE: group_normalization = {group_normalization}')
        if NOTICE:
            print(DEBUG_txt[-1])

        # print(f" exp_name: {exp_name}\n {DataGroup1_name}: {DataGroup1}\n {DataGroup2_name}: {DataGroup2}\n use_LFQ: {use_LFQ}\n useUniProt: {useUniProt}")

# Vnos molekulskih mass
class PeptideVisualizerSetMw(tk.Tk):
    def __init__(self, *args, **kwargs):
        global N_max
        protein_length_example = 255
        tk.Tk.__init__(self, *args, **kwargs)
        
        screen_w = self.winfo_screenwidth()
        screen_h = self.winfo_screenheight()
        screen_dpi = self.winfo_fpixels('1i')
        OS_scale = screen_dpi / 96 # 96 is a scale factor for 100%
        scale_f = screen_dpi / DEV_DPI
        screen_w_eff = round(screen_w / OS_scale)
        screen_h_eff = round(screen_h / OS_scale)
        wh = round(900 * OS_scale)
        ww = round(550 * OS_scale)
        
        DEBUG_txt.append(f"NOTICE: resolution={screen_w}x{screen_h} dpi={screen_dpi} OS_scale={OS_scale} ({screen_w_eff}x{screen_h_eff})")
        DEBUG_txt.append(f"NOTICE: scale_f={scale_f} => {wh}x{ww}")
        if NOTICE:
            print(DEBUG_txt[-2]+"\n"+DEBUG_txt[-1])
        
        self.tk.call('tk', 'scaling', scale_f)
        self.title('Set Y-ax labels')
        self.geometry(f"{wh}x{ww}")
        self.resizable()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [2, 1.2]})
        fig.suptitle("Example peptograph")
        ax1.set_xlim(1, protein_length_example)

        if useUniProt_features:
            ax1.set_ylim(N_max+0.5+(uniProtFeature_coeficient*(N_max+1)), 0.5) # če features, potem k*(N_max+1) večji graf
            ax1.axhline(y=N_max+0.5, color='#000000', linestyle='-', linewidth=0.4)

            #legendbar_height =  uniProtFeature_coeficient*(N_max+1) / len(useUniProt_features)
            #for i in range(len(useUniProt_features)):
            #    for feature in FASTA_features[protein][UniProt_Feature_types[useUniProt_features[i]]]:
            #        ax1.add_patch(Rectangle((feature[1], N_max+0.5+(i*legendbar_height)+(0.1*legendbar_height)), feature[2], 0.8*legendbar_height, alpha=UniProt_Features_alpha, linewidth=0.4, edgecolor='k', facecolor=str(feature[4])))
        else:
            ax1.set_ylim(N_max+0.5, 0.5)

        y_ticks = []
        for i in range(N_max+1):
            y_ticks.append(i+0.5)
        labels = [] * N_max
        ax1.set_yticks(y_ticks, labels)
        
        ax1.set_ylabel('Molecular weight [kDa]')
        ax1.set_yticks(y_ticks, labels)
        ax1.yaxis.set_label_coords(-0.1, 0.5)
        ax1.set_xlabel('Amino acid position', labelpad=0)
        ax1.set_xticks([1, protein_length_example])

        secAx1 = ax1.twiny()
        secAx1.set_xlim(1, protein_length_example)
        secAx1.set_xticks([1, protein_length_example], ["N-term", "C-term"])

        ax2.set_xlim(0, 1.01)
        if useUniProt_features:
            ax2.set_ylim(N_max+0.5+(uniProtFeature_coeficient*(N_max+1)), 0.5)
            ax2.axhline(y=N_max+0.5, color='#000000', linestyle='-', linewidth=0.4)
        else:
            ax2.set_ylim(N_max+0.5, 0.5)

        ax2.set_xlabel('Slice intensity', labelpad=10)
        ax2.tick_params(bottom=False)
        ax2.set_yticks(list(range(1, N_max+1)))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        
        #horizontalne črte
        if horizontal_lines:
            for i in range(N_max):
                ax1.axhline(y=i+0.5, color='#000000', linestyle='-', linewidth=0.1)
                ax2.axhline(y=i+0.5, color='#000000', linestyle='-', linewidth=0.1)
        
        y_offest = 0.1
        plt.subplots_adjust(left = y_offest, bottom = 0.11, right = 0.95, top = 0.87, wspace=0)
        
        canvas = FigureCanvasTkAgg(fig, master = self)  
        canvas.draw()
        canvas.get_tk_widget().place(x=0, y=0)
        
        self.setMWlabels = list()
        self.setMWtexts = list()
        spacing_factor = 355 / N_max
        
        for i in range(N_max-1):
            self.setMWtext = tk.Entry(width=5)
            self.setMWtext.place(x=round(50*OS_scale), y=round((55+((i+1)*spacing_factor))*OS_scale))
            #self.setMWtext.place(x=50, y=(90+(spacing_factor/2)+(spacing_factor*i)))
            self.setMWtexts.append(self.setMWtext)
        
        #presets
        MW_options = []
        for preset_name in MW_PRESET.keys():
            if N_max in MW_PRESET[preset_name]:
                MW_options.append(preset_name)
        if len(MW_options) != 0:
            self.preset_name = tk.StringVar(value="Gel presets")
            self.preset_menu = tk.OptionMenu(self, self.preset_name, *MW_options, command=self.set_preset).place(x=round(100*OS_scale), y=round(507*OS_scale))
        
        self.setYax_btn = tk.Button(text="Submit", command=self.setYax).place(x=round(420*OS_scale), y=round(510*OS_scale))
    
    def set_preset(self, arrg): # od kje dobi drugi argument??
        del arrg
        global N_max
        #print(MW_PRESET[self.preset_name.get()][N_max])
        for i in range(N_max-1):
            self.setMWtexts[i].delete(0, tk.END)
            self.setMWtexts[i].insert(0, MW_PRESET[self.preset_name.get()][N_max][i])
        
    def setYax(self):
        global N_max
        global setMW_Labels
        for i in range(N_max-1):
            setMW_Labels.append(self.setMWtexts[i].get())
        self.quit()
        self.destroy()

suffixes = ['LFQ intensity ', 'Intensity ']
patterns = [re.compile(suffix) for suffix in suffixes]

def remove_obvious(s: str) -> str:
    for pattern in patterns:
        s = pattern.sub("", s)
    return s

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return False

def filelinecount(diz_file):
    count = 0
    thefile = open(diz_file)
    while 1:
        buffer = thefile.read(65536)
        if not buffer: break
        count += buffer.count('\n')
    thefile.close()
    del thefile
    return count

def FastaHeader2NameID(tmp):
    tmp = tmp.strip()
    #poberem sam do OS
    if tmp.find('OS=') != -1:
        tmp = tmp[:tmp.find('OS=')]
    
    tmp = tmp.split('|')
    if len(tmp) > 1:
        if tmp[2].find(" ") > 0:
            FastaHeader2NameID_name = tmp[2][tmp[2].find(" "):]
            FastaHeader2NameID_id = tmp[2][:tmp[2].find(" ")]
        return FastaHeader2NameID_name, FastaHeader2NameID_id
    else:
        
        return tmp[0],False
        
# UniProt magic
def UniProt_FASTA(identifier):
    #HTTP_request_fasta = requests.get(f'https://www.uniprot.org/uniprot/{identifier}.fasta')
    HTTP_request_fasta = requests.get(f'https://rest.uniprot.org/uniprotkb/{identifier}.fasta')

    if(HTTP_request_fasta.status_code == 200):
        diz_OUT = ''
        seq_start = False
        for line in HTTP_request_fasta.text.splitlines():
            if seq_start:
                diz_OUT += line.strip()
            else:
                diz_header = line.strip()
                seq_start = True
    else:
        print(f"FATAL ERROR! UniProt responded with code {HTTP_request_fasta.status_code}! Unable to retrive FASTA ({identifier})!")
        exit(2) # a to sploh vrže vn???
    del HTTP_request_fasta
    return diz_header, diz_OUT

# USAGE: FASTA_features[identifier] = UniProt_Features(identifier)
def UniProt_Features(identifier):
    #HTTP_request_fasta = requests.get(f'https://www.uniprot.org/uniprot/{identifier}.gff')
    HTTP_request_fasta = requests.get(f'https://rest.uniprot.org/uniprotkb/{identifier}.gff')
    
    if(HTTP_request_fasta.status_code == 200):
        this_feaures = {}
        this_feaures[UniProt_Feature_types[0]] = []
        this_feaures[UniProt_Feature_types[1]] = []
        this_feaures[UniProt_Feature_types[2]] = []
        for line in HTTP_request_fasta.text.splitlines():
            feature_fields = re.split(r'\t+', line.rstrip('\t'))
            if len(feature_fields) >= 8: # tako zaznam da ni header ali kaj podobnega
                #filter za feature
                feature_length = int(feature_fields[4]) - int(feature_fields[3])
                if feature_length > 1 and feature_length < len(FASTA[identifier]) - 1:
                    note_result = find_between(feature_fields[8], 'Note=', ';') # a ima feature kako opombo?
                    if any(feature_fields[2] in x for x in LEGEND):
                        x = [x for x in LEGEND if feature_fields[2] in x][0]
                        color = x[1]
                    else:
                        color = LEGEND[0][1]
                        
                    if(feature_fields[2] in UniProt_Processing_Features): # Molecular processing
                        this_feaures[UniProt_Feature_types[0]].append([feature_fields[2], int(feature_fields[3]), feature_length, note_result, color])
                    if(feature_fields[2] in UniProt_Region_Features): # Topology / region
                        this_feaures[UniProt_Feature_types[1]].append([feature_fields[2], int(feature_fields[3]), feature_length, note_result, color])
                    if(feature_fields[2] in UniProt_2ndary_Features): # Secondary structure
                        this_feaures[UniProt_Feature_types[2]].append([feature_fields[2], int(feature_fields[3]), feature_length, note_result, color])
            #print(feature_fields)
    else:
        print(f"ERROR! UniProt responded with code {HTTP_request_fasta.status_code}! Unable to retrive features! Terminating ...")
        pause()
        exit(2)
    del HTTP_request_fasta
    return(this_feaures)    

def write_debug():
    print("Writing debug ...")
    
    debug_file = "debug_"+datetime.now().strftime("%Y-%m-%d_%H-%M")+".txt"
    
    if os.path.exists(debug_file):
        os.remove(debug_file)
    file = open(debug_file, "w")
    # DEBUG_txt
    for debug_line in DEBUG_txt:
        file.write(debug_line+"\n\n")
    # replicates
    file.write("replicates: ")
    try:
        replicates
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(str(replicates)+"\n")
    # N_max
    file.write("N_max: ")
    try:
        N_max
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(str(N_max)+"\n")
    # peptide_threshold
    file.write("peptide_threshold: ")
    try:
        peptide_threshold
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(str(peptide_threshold*100)+" %\n")
    # data_titles
    file.write("data_titles: ")
    try:
        data_titles
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(data_titles)+"\n")
    # data_indexes
    file.write("data_indexes: ")
    try:
        data_indexes
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(data_indexes)+"\n")
    # DataGroup1
    file.write("DataGroup1: ")
    try:
        DataGroup1
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(DataGroup1)+"\n")
    # DataGroup2
    file.write("DataGroup2: ")
    try:
        DataGroup2
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(DataGroup2)+"\n")
    # SUMMARY
    file.write("SUMMARY: ")
    try:
        SUMMARY
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(SUMMARY)+"\n")
    # DB
    file.write("DB: ")
    try:
        DB
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(DB)+"\n")
    # K_factor
    file.write("K_factor: ")
    try:
        K_factor
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(K_factor)+"\n")
    # NOZMALIZATION_FACTORS
    file.write("NOZMALIZATION_FACTORS: ")
    try:
        NOZMALIZATION_FACTORS
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(NOZMALIZATION_FACTORS)+"\n")
    # SLICE_quantitative
    file.write("SLICE_quantitative: ")
    try:
        SLICE_quantitative
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(SLICE_quantitative)+"\n")
    # SLICE_quantitative_err
    file.write("SLICE_quantitative_err: ")
    try:
        SLICE_quantitative_err
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(SLICE_quantitative_err)+"\n")
    # PEPTIDES_qualitative
    file.write("PEPTIDES_qualitative: ")
    try:
        PEPTIDES_qualitative
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(PEPTIDES_qualitative)+"\n")
    # PEPTIDES_quantitative
    file.write("PEPTIDES_quantitative: ")
    try:
        PEPTIDES_quantitative
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(PEPTIDES_quantitative)+"\n")
    # PEPTIDES_quantitative_avg
    file.write("PEPTIDES_quantitative_avg: ")
    try:
        PEPTIDES_quantitative_avg
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(PEPTIDES_quantitative_avg)+"\n")
    # PEPTIDES_quantitative_err
    file.write("PEPTIDES_quantitative_err: ")
    try:
        PEPTIDES_quantitative_err
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(PEPTIDES_quantitative_err)+"\n")
    # PEPTIDES_quantitative_avg_norm
    file.write("PEPTIDES_quantitative_avg_norm: ")
    try:
        PEPTIDES_quantitative_avg_norm
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(PEPTIDES_quantitative_avg_norm)+"\n")
    # FASTA_name
    file.write("FASTA_name: ")
    try:
        FASTA_name
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(FASTA_name)+"\n")
    # FASTA
    file.write("FASTA: ")
    try:
        FASTA
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(FASTA)+"\n")
    # FASTA_features
    file.write("FASTA_features: ")
    try:
        FASTA_features
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(FASTA_features)+"\n")

    file.close()

###########
#
# ONLINE CHECK + UPDATE MODUL
#
###########

for i in range(len(online_srv)):
            
    ping_arg = "-n 1" if OS=="Windows" else "-c 1"
    ping_cmd = "ping "+ping_arg+" "+online_srv[i] if OS=="MacOS" else "ping "+ping_arg+" -w 2 "+online_srv[i]

    DEBUG_txt.append(f'NOTICE: Checking connectivity using command "{ping_cmd}"')
    if NOTICE:
        print(DEBUG_txt[-1])

    status,result = subprocess.getstatusoutput(ping_cmd)
    if status == 0:
        ONLINE = True
        break
    else:
        ONLINE = False
        DEBUG_txt.append(f"NOTICE: No respons from {online_srv[i]}")
        if NOTICE:
            print(DEBUG_txt[-1])

if not ONLINE:
    DEBUG_txt.append(f"WARNING! {PROGRAM_NAME} is unable to access internet! Skipping updates check. UniProt features disabled.")
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")

# ONLINE
DEBUG_txt.append(f'NOTICE: ONLINE = {ONLINE}')
if NOTICE:
    print(DEBUG_txt[-1])
# SELFUPDATE
DEBUG_txt.append(f'NOTICE: SELFUPDATE = {SELFUPDATE}')
if NOTICE:
    print(DEBUG_txt[-1])

UPDATED = False
if SELFUPDATE and ONLINE:
    # DEJANSKI APDEJT
    DEBUG_txt.append(f"NOTICE: Checking for update @ \"{update_url}\"")
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        HTTP_UPDATE = requests.get(update_url)
        if(HTTP_UPDATE.status_code == 200):
            online_mtime = int(calendar.timegm(parsedate(HTTP_UPDATE.headers['Last-Modified']).timetuple())) # calendar.timegm() namesto time.mktime() ker je mtime v HTTP headerju podan v UTC, ne v local
            local_mtime = int(os.path.getmtime(__file__))
            DEBUG_txt.append(f"NOTICE: local mtime: {datetime.fromtimestamp(local_mtime).strftime('%Y-%m-%d %H:%M:%S')}, online mtime: {datetime.fromtimestamp(online_mtime).strftime('%Y-%m-%d %H:%M:%S')}")
            if NOTICE:
                print(DEBUG_txt[-1])
            
            if online_mtime > local_mtime:
                DEBUG_txt.append(f"WARNING! New version of {PROGRAM_NAME} detected!")
                print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                if(query_yes_no(f'Do you want to update {PROGRAM_NAME}?')):
                    
                    try:
                        shutil.copy2(__file__, __file__+'.v'+VERSION+'.old')
                        with open(__file__, 'wb') as f:
                            f.write(HTTP_UPDATE.content)

                    except:
                        DEBUG_txt.append(f"ERROR! Update failed!")
                        print(f"{bcolors.FAIL}{DEBUG_txt[-1]}{bcolors.ENDC}")
                        pause()
                        exit(2)
                
                    UPDATED = True
                else:
                    DEBUG_txt.append(f"NOTICE: Update skipped.")
                    if NOTICE:
                        print(DEBUG_txt[-1])

            else:
                DEBUG_txt.append(f"NOTICE: local version is up to date.")
                if NOTICE:
                    print(DEBUG_txt[-1])
        else:
            DEBUG_txt.append(f"NOTICE: {PROGRAM_NAME} is unable to access update server (SRV CODE: {HTTP_UPDATE.status_code})!")
            if NOTICE:
                print(DEBUG_txt[-1])
    except:
        DEBUG_txt.append(f"NOTICE: {PROGRAM_NAME} is unable to access update server (UNKNOWN REASON)!")
        if NOTICE:
            print(DEBUG_txt[-1])

if UPDATED:
    DEBUG_txt.append(f"SUCCESS! {PROGRAM_NAME} has been successfully updated! Please restart script.")
    print(f"{bcolors.OKGREEN}{DEBUG_txt[-1]}{bcolors.ENDC}")
    pause()
    exit()

##############################################################################
#
#                       I. Uporabniški vnosi
#
##############################################################################

# uporabnik izbere fajl peptides.txt

DEBUG_txt.append('Waiting for peptides.txt ... ')
if GUI:
    print(DEBUG_txt[-1])
    peptides_file = askopenfilename(title = "Select peptides file",filetypes=[("peptides.txt", ".txt .csv")])
else:
    print("")
    peptides_file = input("Enter absolute location of peptides.txt: ").strip()

if not peptides_file or not exists(peptides_file):
    DEBUG_txt.append("ERROR! No peptides.txt file provided. Terminating ...")
    print(f"{bcolors.FAIL}{DEBUG_txt[-1]}{bcolors.ENDC}")
    exit(2)
DEBUG_txt.append(f'NOTICE: Using file "{peptides_file}".')
if NOTICE:
    print(DEBUG_txt[-1])

# v1.2 pridobivanje data_titles iz prve vrstice preden ponudim možnost LFQ
DEBUG_txt.append('Reading first line of peptides.txt ...')
print(DEBUG_txt[-1])

existsLFQ = True
with open(peptides_file) as raw_peptides_file:
    csv_reader = csv.reader(raw_peptides_file, delimiter="\t")
    for row in csv_reader:
        ###
        ##### v prvi vrstici ugotovim kateri stolpci me zanimajo + uporabnik pove kateri zanimajo njega #####
        ###
        sequence_index = row.index('Sequence')
        protein_index = row.index('Leading razor protein')
        start_position_index = row.index('Start position')
        length_index = row.index('Length')
        contaminant_index = row.index('Potential contaminant')
        reverse_index = row.index('Reverse')
        #LFQ intensity/Intensity?
        
        data_titles_LFQ = [i for i in row if i.startswith('LFQ intensity ')]
        if not data_titles_LFQ:
            existsLFQ = False
            DEBUG_txt.append(f"NOTICE: LFQ data not found ... ")
            if NOTICE:
                print(f"\n{DEBUG_txt[-1]}")
        data_titles_INT = [i for i in row if i.startswith('Intensity ')]
            
        if not data_titles_INT:
            print(f"{bcolors.FAIL}ERROR! No data detected!{bcolors.ENDC}")
            pause()
            exit()
        
        if existsLFQ:
            data_indexes_LFQ = []
            for data_title in data_titles_LFQ:
                data_indexes_LFQ.append(row.index(data_title))
            
        data_indexes_INT = []
        for data_title in data_titles_INT:
            data_indexes_INT.append(row.index(data_title))
            
        #odstranim nepotrebne dodatke imen
        if existsLFQ:
            data_titles_LFQ = [w.replace('LFQ intensity ', '') for w in data_titles_LFQ]
        data_titles_INT = [w.replace('Intensity ', '') for w in data_titles_INT]
        break # samo prva vrstica me zanima

# GUI jupiii
if GUI:
    #while not DataGroup1:
    App = PeptideVisualizerGUI()
    App.mainloop()
    if not DataGroup1:
        print(f"{bcolors.FAIL}ERROR! Experiment settings not set! Terminating ...{bcolors.ENDC}")
        pause()
        exit()
else:
    # LFQ / intensity?
    if existsLFQ:
        print("")
        use_LFQ = query_yes_no('Do you want to use LFQ intensities?', 'no')
        DEBUG_txt.append(f'NOTICE: use_LFQ = {use_LFQ}')
        if NOTICE:
            print(DEBUG_txt[-1])
    else:
        use_LFQ = False

    if use_LFQ:
        data_indexes = data_indexes_LFQ
        data_titles = data_titles_LFQ
    else:
        data_indexes = data_indexes_INT
        data_titles = data_titles_INT

    # UniProt?

    print("")
    if not ONLINE:
        useUniProt = False
    else:
        useUniProt = query_yes_no("Use UniProt online database? This will alow feature visualization.")
    DEBUG_txt.append(f'NOTICE: useUniProt = {useUniProt}')
    if NOTICE:
        print(DEBUG_txt[-1])

    if useUniProt:
        useUniProt_features = userUniProtFeatures()
    else:
        useUniProt_features = False
        DEBUG_txt.append(f'NOTICE: useUniProt_features = {useUniProt_features}')
        if NOTICE:
            print(DEBUG_txt[-1])

    ##############################################################################
    #                Uporabnik nastavi skupine, način normalizacije in ostalo
    ##############################################################################

    DataGroup1, DataGroup2 = userGrouping(data_titles)
    N_max = len(DataGroup1)

    # Replicates
    print("")
    try:
        replicates = int(input("Number of replicates used [default: 1]: ").replace(" ", ""))
    except ValueError:
        print(f"{bcolors.WARNING}WARNING! Number of replicates could not be interpreted. Using default value 1{bcolors.ENDC}")
        replicates = 1
    if replicates < 1 or replicates > replicates_MAX:
        print(f"{bcolors.WARNING}WARNING! Number of replicates makes no sense ({replicates}). Using default value 1{bcolors.ENDC}")
        replicates = 1
    if N_max % replicates != 0:
        print(f"{bcolors.FAIL}ERROR! Number of replicates ({replicates}) makes no sense considering group size ({N_max}). {N_max} is not divisible by {replicates}! Terminating ...{bcolors.ENDC}")
        pause()
        exit(2)
    
    N_max = int(N_max / replicates) # nov maximum
    DEBUG_txt.append(f"NOTICE: replicates = {replicates}\t => \tN slices: {N_max}")
    if NOTICE:
        print(DEBUG_txt[-1])
    
    # Stilistični vprašanji
    if replicates > 1:
        print("")
        error_bar = query_yes_no("Do you want to display error on slice intensity bars?")
    else:
        error_bar = False
    print("")
    horizontal_lines = query_yes_no("Do you want to display horizontal lines between slices?")

    # Ročno poimenovanje Y osi
    print("")
    setMW_Labels = []
    setMW = query_yes_no('Do you want to label Y-axis with molar mass?', 'no')
    if setMW:
        for n in range(N_max-1):
            try:
                tmp_in = input(f"Molecular weight [kDa] between slice {n+1} and slice {n+2}: ").strip().lower().replace(" ", "").replace("kda", "")
                if not tmp_in:
                    setMW_Labels.append("")
                else:
                    setMW_Labels.append(int(tmp_in))
                del tmp_in
            except ValueError:
                DEBUG_txt.append(f"WARNING! Molecular weight [kDa] between slice {n+1} and slice {n+2} (n={n}) could not be interpreted. Turning custom Y-axis label off.")
                print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                setMW = False
                break;
            

    DEBUG_txt.append(f"NOTICE: setMW = {setMW}. setMW_Labels = {setMW_Labels}")
    if NOTICE:
        print(DEBUG_txt[-1])
    
    # Normalization method
    if DataGroup2:
        print("")
        group_normalization = query_yes_no(f'Do you want to normalize peptide intensities in each group? {bcolors.WARNING}USING THIS IS NOT RECOMMENDED!{bcolors.ENDC}', 'no')

    if group_normalization:
        DEBUG_txt.append('WARNING! Group normalization is used insted of overall normalization!')
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")

    DEBUG_txt.append(f'NOTICE: group_normalization = {group_normalization}')
    if NOTICE:
        print(DEBUG_txt[-1])
    # peptide_threshold
    print("")
    try:
        peptide_threshold = float(input("Relative peptide detection threshold [default: {}%]: ".format(default_peptide_threshold)).replace(",", ".").replace(" ", "").replace("%", ""))
    except ValueError:
        DEBUG_txt.append(f"WARNING! Relative peptide detection threshold could not be interpreted. Using default value {default_peptide_threshold}%")
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
        peptide_threshold = default_peptide_threshold
    if peptide_threshold < 0 or peptide_threshold > 100:
        DEBUG_txt.append(f"WARNING! Relative peptide detection threshold makes no sense ({peptide_threshold}%). Using default value {default_peptide_threshold}%")
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
        peptide_threshold = default_peptide_threshold

    DEBUG_txt.append(f'NOTICE: peptide_threshold = {peptide_threshold}')
    if NOTICE:
        print(DEBUG_txt[-1])

    # Is SVG fine?
    print("")
    SVG_OUTPUT = query_yes_no("Do you want to generate vector (.svg) output?")
    
    DEBUG_txt.append(f'NOTICE: SVG_OUTPUT = {SVG_OUTPUT}')
    if NOTICE:
        print(DEBUG_txt[-1])

    # Ime eksperimenta, ko že veš kakšne nastavitve si nastavil
    print("")
    usr_name = slugify(input("Enter experiment name: "))

if not useUniProt:
    # OFFLINE MODE
    DEBUG_txt.append(f'Waiting for FASTA file due to OFFLINE mode ...')
    if GUI:
        print(DEBUG_txt[-1])
        fasta_file = askopenfilename(title = "Select fasta file",filetypes=[("FASTA file", ".fasta")])
    else:
        print("")
        fasta_file = input("Enter absolute location of FASTA file: ").strip()

    if not fasta_file or not exists(fasta_file):
        print(f"{bcolors.FAIL}ERROR! UniProt online database dissabled and no FASTA file provided. Terminating ...{bcolors.ENDC}")
        pause()
        exit(2)
    DEBUG_txt.append(f'NOTICE: Using offline FASTA file "{fasta_file}"')
    if NOTICE:
        print(DEBUG_txt[-1])

# ALL INPUT DONE
DEBUG_txt.append('NOTICE: ALL USER INPUTS SUCCESSFUL!')
if NOTICE:
    print(DEBUG_txt[-1])

if not usr_name:
    usr_name = PROGRAM_NAME

now = datetime.now()
exp_name = now.strftime("%Y-%m-%d_%H-%M_") + usr_name
del now
work_dir = os.path.dirname(peptides_file) + '/' + exp_name + '/'
DEBUG_txt.append(f'NOTICE: Experiment name: {usr_name}. Assuming working directory: "{work_dir}"')
if NOTICE:
    print(DEBUG_txt[-1])

# K factor threshold (in # of slices)
K_factor_threshold = int(math.ceil(N_max*0.05)) # 1 do 20 slicov, drugače 2
DEBUG_txt.append(f'NOTICE: K_factor_threshold = {K_factor_threshold}')
if NOTICE:
    print(DEBUG_txt[-1])

if not DataGroup2:
    group_normalization = False

peptide_threshold *= 0.01

if SVG_OUTPUT:
    OUTPUT_EXTENTION = '.svg'
else:
    OUTPUT_EXTENTION = '.png'

##############################################################################
#               peptides.txt handling
##############################################################################
start_time = time.time()
FASTA = {}
FASTA_name = {}
FASTA_features = {}

DEBUG_txt.append('Reading and filtering peptides.txt ...')
print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")
# v DB bo ključ protein in peptidi seznam
DB = {}
# tukaj bo ključ peptid in intenzitete seznam
PEPTIDES_quantitative = {}
# tukaj bo ključ peptid in bodo napisani peptidi od kje do kje
PEPTIDES_qualitative = {}

###
##### Pridobivanje podatkov #####
###
header_conflict = False
with open(peptides_file) as raw_peptides_file:
    csv_reader = csv.reader(raw_peptides_file, delimiter="\t")
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            tmp_protein = row[protein_index]
            if row[reverse_index] != '+' and row[contaminant_index] != '+' and tmp_protein[0:5] != 'REV__': #Tretji stavek je fail-saffe, ker MQ včasih ne doda + v reverse
                if useUniProt:
                    tmp_protein_identifier = find_between(tmp_protein, '|', '|')
                    if not tmp_protein_identifier:
                        if not header_conflict:
                            DEBUG_txt.append(f"WARNING: peptides.txt does not use a valid UniProt header ({tmp_protein})! Assuming UniProt identifier ...")
                            print(bcolors.WARNING+DEBUG_txt[-1]+bcolors.ENDC)
                            header_conflict=True
                        tmp_protein_identifier = tmp_protein
                
                tmp_peptide = row[sequence_index]
                tmp_start = int(row[start_position_index]) - 1
                tmp_length = int(row[length_index])
                PEPTIDES_qualitative[tmp_peptide] = [tmp_start, tmp_length]
                
                tmp_data = []
                for index in data_indexes:
                    tmp_data.append(float(row[index]))
                
                if useUniProt:
                    DB.setdefault(tmp_protein_identifier, []).append(tmp_peptide)
                else:
                    DB.setdefault(tmp_protein, []).append(tmp_peptide)
                PEPTIDES_quantitative[tmp_peptide] = tmp_data
                #just in case
                if useUniProt:
                    del tmp_protein_identifier
                del tmp_protein
                del tmp_peptide
                del tmp_start
                del tmp_length
                del tmp_data
            elif tmp_protein[0:5] == 'REV__' and row[reverse_index] != '+':
                DEBUG_txt.append(f'{bcolors.WARNING}WARNING: MaxQuant did not assign {tmp_protein} as reverse!{bcolors.ENDC}')
                print(DEBUG_txt[-1])
        line_count += 1
    
read_time = time.time() - start_time
print(f"{bcolors.OKBLUE}Read {line_count} peptides in {round(read_time*1000)} ms{bcolors.ENDC}")

#prepare list for processing
multi_IDs = []
for id in DB:
    multi_IDs.append(id) 

#print(multi_IDs)
#pause()

DEBUG_txt.append(f'All the required data is stored in the memory. Starting peptide visualization of {len(DB.keys())} proteins ...')
if NOTICE:
    print(DEBUG_txt[-1])

##############################################################################
#               FASTA handling
##############################################################################
continue_time = time.time()
if useUniProt:
    DEBUG_txt.append(f" ( 1/3 ) Retriving FASTA from UniProt [{continue_time}] ...")
    print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")

    if Multiproc:
        # Multiprocessing
        def getUniProt(id):
            header, FASTA[id] = UniProt_FASTA(id)
            FASTA_name[id], header_id = FastaHeader2NameID(header)
            if useUniProt_features:
                FASTA_features[id] = UniProt_Features(id)
            return id
        
        with alive_bar(len(multi_IDs), dual_line=True, theme='classic') as bar:
            results = ThreadPool(N_download_threads).imap_unordered(getUniProt, multi_IDs)
            for r in results:
                bar.text = f'Retrived UniProt data for {r}'
                bar()
    else:
        # SINGLE THREAD
        with alive_bar(len(DB.keys()), dual_line=True, theme='classic') as bar:
            for identifier in DB:
                bar.text = f'Retriving UniProt data for identifier {identifier}'
            
                header, FASTA[identifier] = UniProt_FASTA(identifier)
                FASTA_name[identifier], header_id = FastaHeader2NameID(header)
                if useUniProt_features:
                    FASTA_features[identifier] = UniProt_Features(identifier)
                bar()
                
else:
    DEBUG_txt.append(f" ( 1/3 ) Reading FASTA [{continue_time}] ...")
    print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")
    with open(fasta_file) as raw_fasta_file:
        for line in raw_fasta_file:
            #print(line)
            if line[0] == '>':
                sp=line[1:].split(' ')
                name = sp[0]
                FASTA[name] = ''
                ## ERR TEST REQUIRED
                FASTA_name[name], header_id = FastaHeader2NameID(line)
            else:
                FASTA[name] += line.strip()


##############################################################################
#               Statistična obdelava rezultatov
##############################################################################

DEBUG_txt.append(f" ( 2/3 ) Statistically analyzing quantities and normalizing them [{time.time()}] ... ")
print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")

NOZMALIZATION_FACTORS = {}
PEPTIDES_quantitative_avg = {}
PEPTIDES_quantitative_err = {}
PEPTIDES_quantitative_avg_norm = {}
SLICE_quantitative = {}
SLICE_quantitative_err = {}
Total_intensity_1 = {}
Total_intensity_2 = {}
K_factor = {}
SUMMARY = {}

with alive_bar(len(DB.keys()), dual_line=True, theme='classic') as bar:
    for this_protein, these_peptides in DB.items():
        bar.text = f'Statistically analyzing protein {this_protein}'
        #FASTA check
        try: # to ne dela, ker se sesuje. morda bi lahko prečekiral pri fasti, če je kompletna?
            FASTA[this_protein]
        except NameError:
            print(f'{bcolors.FAIL}ERROR! FASTA entry for "{this_protein}" does not exist! Check FASTA file. Terminating ...{bcolors.ENDC}')
            pause()
            exit(2)
        ###
        # SLICE_quantitative in SLICE_quantitative_err, dokler še imam intenzitete absolutne.
        ###

        SLICE_quantitative[this_protein] = {}
        SLICE_quantitative_err[this_protein] = {}
        slice = 1
    
        for i in range(N_max):
            ID = DataGroup1[i]
            replicate_intensity = 0
            work_list = []

            for j in range(replicates):
                rep_ID = DataGroup1[i+(j*N_max)]
                for this_peptide in these_peptides:
                    replicate_intensity += PEPTIDES_quantitative[this_peptide][rep_ID]
                    
                work_list.append(replicate_intensity)
                    
            if replicates > 1:
                SLICE_quantitative[this_protein][slice] = statistics.mean(work_list)
                SLICE_quantitative_err[this_protein][slice] = statistics.stdev(work_list)
            else:
                SLICE_quantitative[this_protein][slice] = work_list[0]
                SLICE_quantitative_err[this_protein][slice] = 0
            del work_list
            slice += 1
            
        if DataGroup2:
            for i in range(N_max):
                ID = DataGroup2[i]
                replicate_intensity = 0
                work_list = []
            
                for j in range(replicates):
                    rep_ID = DataGroup2[i+(j*N_max)]
                    for this_peptide in these_peptides:
                        replicate_intensity += PEPTIDES_quantitative[this_peptide][rep_ID]
                    
                    work_list.append(replicate_intensity)
                        
                if replicates > 1:
                    SLICE_quantitative[this_protein][slice] = statistics.mean(work_list)
                    SLICE_quantitative_err[this_protein][slice] = statistics.stdev(work_list)
                else:
                    SLICE_quantitative[this_protein][slice] = work_list[0]
                    SLICE_quantitative_err[this_protein][slice] = 0
                del work_list
                slice += 1
        
        # skupna intenziteta
        Total_intensity_1[this_protein] = Total_intensity_2[this_protein] = 0
        for i in range(N_max):
            Total_intensity_1[this_protein] += SLICE_quantitative[this_protein][i+1]
            if DataGroup2:
                Total_intensity_2[this_protein] += SLICE_quantitative[this_protein][N_max+i+1]
        
        tmp_norm = max(Total_intensity_1[this_protein], Total_intensity_2[this_protein])
        if tmp_norm != 0:
            Total_intensity_1[this_protein] /= tmp_norm
            Total_intensity_2[this_protein] /= tmp_norm
        del tmp_norm

        # normalizacija rezin
        slice_normalization_factor = max(SLICE_quantitative[this_protein].values())
        if slice_normalization_factor > 0:
            for n in range(len(SLICE_quantitative[this_protein].keys())):
                SLICE_quantitative[this_protein][n+1] = SLICE_quantitative[this_protein][n+1] / slice_normalization_factor
                SLICE_quantitative_err[this_protein][n+1] = SLICE_quantitative_err[this_protein][n+1] / slice_normalization_factor
        del slice_normalization_factor
        
        ###
        # Izračun k-faktorja za sortiranje. Kasneje je pomnožen še z pokritostjo.
        ###
        if DataGroup2:
            K_factor[this_protein] = 0
            for i in range(N_max):
                for j in range(N_max):
                    if i < j: # da se primerjava ne ponavlja
                        distance = abs(i-j)-K_factor_threshold
                        delta_i = abs( SLICE_quantitative[this_protein][N_max+i+1] - SLICE_quantitative[this_protein][i+1] )
                        delta_j = abs( SLICE_quantitative[this_protein][N_max+j+1] - SLICE_quantitative[this_protein][j+1] )
                        if distance > 0: # če je razlika med i in j več kot threshold
                            K_factor[this_protein] += 100 * distance * (delta_i * delta_j) + delta_i + delta_j
        
        
        ###
        # transformiram PEPTIDES_quantitative v PEPTIDES_quantitative_avg in PEPTIDES_quantitative_err pa pol uporabljam ta dva kasneje bojo normalizirani
        ###
        for this_peptide in these_peptides:
            if replicates > 1:
                PEPTIDES_quantitative_avg[this_peptide] = {}
                PEPTIDES_quantitative_err[this_peptide] = {}
                PEPTIDES_quantitative_avg_norm[this_peptide] = {}
        
                avg_list_1 = []
                err_list_1 = []
                avg_list_2 = []
                err_list_2 = []
            
                for i in range(N_max):
                    #DataGroup1
                    ID = DataGroup1[i]
                    work_list = []
                    for j in range(replicates):
                        rep_ID = DataGroup1[i+(j*N_max)]
                        #work_list.append(PEPTIDES_quantitative[this_peptide][ID+(j*N_max)]) FAIL
                        work_list.append(PEPTIDES_quantitative[this_peptide][rep_ID])
                        
                    avg_list_1.append(statistics.mean(work_list))
                    if avg_list_1[i] > 0: # da ne pride do deljenja z 0
                        err_list_1.append(statistics.stdev(work_list))
                    else:
                        err_list_1.append(0)
                    del(work_list)
                    #DataGroup2
                    if DataGroup2:
                        ID = DataGroup2[i]
                        work_list = []
                        for j in range(replicates):
                            rep_ID = DataGroup2[i+(j*N_max)]
                            work_list.append(PEPTIDES_quantitative[this_peptide][rep_ID])
                        avg_list_2.append(statistics.mean(work_list))
                        if avg_list_2[i] > 0: # da ne pride do deljenja z 0
                            err_list_2.append(statistics.stdev(work_list))
                        else:
                            err_list_2.append(0)
                        del(work_list)
                if DataGroup2:
                    PEPTIDES_quantitative_avg[this_peptide] = avg_list_1+avg_list_2
                    PEPTIDES_quantitative_err[this_peptide] = err_list_1+err_list_2
                else:
                    PEPTIDES_quantitative_avg[this_peptide] = avg_list_1
                    PEPTIDES_quantitative_err[this_peptide] = err_list_1
                del avg_list_1
                del err_list_1
                del avg_list_2
                del err_list_2
            else: # če ni replikatov je avg samo prepisana vrednost
                work_list = []
                for i in range(N_max):
                    ID = DataGroup1[i]
                    work_list.append(PEPTIDES_quantitative[this_peptide][ID])
                if DataGroup2:
                    for i in range(N_max):
                        ID = DataGroup2[i]
                        work_list.append(PEPTIDES_quantitative[this_peptide][ID])
                        
                PEPTIDES_quantitative_avg[this_peptide] = work_list
                del work_list
                if DataGroup2:
                    PEPTIDES_quantitative_err[this_peptide] = [0] * (2*N_max)
                else:
                    PEPTIDES_quantitative_err[this_peptide] = [0] * N_max
        
        ###
        # Normalization of _avg
        ###
        if group_normalization:
            tmp_list=[]
            peptides_count = 0
            for this_peptide in these_peptides:
                for i in range(N_max):
                    tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][i])
                peptides_count+=1
            #SUMMARY
            SUMMARY[this_protein] = [slugify(this_protein.replace("|", "_")),sum(tmp_list),peptides_count]
        
            normalization_factor = max(tmp_list)
            if normalization_factor > 0:
                NOZMALIZATION_FACTORS[this_protein] = normalization_factor
                for this_peptide in these_peptides:
                    work_list = []
                    for i in range(N_max):
                        work_list.append(PEPTIDES_quantitative_avg[this_peptide][i] / normalization_factor)
                    PEPTIDES_quantitative_avg_norm[this_peptide] = work_list
                    del work_list
            else:
                NOZMALIZATION_FACTORS[this_protein] = False
                for this_peptide in these_peptides:
                    PEPTIDES_quantitative_avg_norm[this_peptide] = [0] * N_max
                
            del normalization_factor
            del tmp_list
            #DataGroup2
            tmp_list=[]
            for this_peptide in these_peptides:
                for i in range(N_max):
                    tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][i])
            #SUMMARY
            SUMMARY[this_protein][1] = sum(tmp_list, SUMMARY[this_protein][1])
        
            normalization_factor = max(tmp_list)
            if normalization_factor > 0:
                for this_peptide in these_peptides:
                    work_list = []
                    for i in range(N_max):
                        work_list.append(PEPTIDES_quantitative_avg[this_peptide][i] / normalization_factor)
                    PEPTIDES_quantitative_avg_norm[this_peptide].extend(work_list) # dodajam k seznamu
                    del work_list
            else:
                for this_peptide in these_peptides:
                    PEPTIDES_quantitative_avg_norm[this_peptide].extend([0] * N_max) #dodajam k seznamu
            del normalization_factor
            del tmp_list
        else: # Normalna normalizacija
            tmp_list=[]
            peptides_count = 0
            for this_peptide in these_peptides:
                peptides_count+=1
                
                for i in range(N_max):
                    tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][i])
                if DataGroup2:
                    for i in range(N_max):
                        tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][N_max+i])
            #SUMMARY
            SUMMARY[this_protein] = [slugify(this_protein.replace("|", "_")),sum(tmp_list),peptides_count]
            
            normalization_factor = max(tmp_list)
            if normalization_factor > 0:
                NOZMALIZATION_FACTORS[this_protein] = normalization_factor
                for this_peptide in these_peptides:
                    work_list = []
                    for i in range(N_max):
                        work_list.append(PEPTIDES_quantitative_avg[this_peptide][i] / normalization_factor)
                    if DataGroup2:
                        for i in range(N_max):
                            work_list.append(PEPTIDES_quantitative_avg[this_peptide][N_max+i] / normalization_factor)
                    PEPTIDES_quantitative_avg_norm[this_peptide] = work_list
                    del work_list
            else: #če je normalčizacijski faktor enak 0 pomeni da so vsi elementi enaki 0
                NOZMALIZATION_FACTORS[this_protein] = False
                if DataGroup2:
                    work_list = [0] * (2*N_max)
                else:
                    work_list = [0] * (N_max)
                for this_peptide in these_peptides:
                    PEPTIDES_quantitative_avg_norm[this_peptide] = work_list
                del work_list
            del normalization_factor
            del tmp_list
        
        bar() #progress

###
##### Generiranje .TXT in .SVG rezultatov v odvisnosti od načina (align/protomap) / srce skripte
###
#
#   (c) Unknown author
#
#          |  \ \ | |/ /
#          |  |\ `' ' /
#          |  ;'aorta \      / , pulmonary
#          | ;    _,   |    / / ,  arteries
# superior | |   (  `-.;_,-' '-' ,
#vena cava | `,   `-._       _,-'_
#          |,-`.    `.)    ,<_,-'_, pulmonary
#         ,'    `.   /   ,'  `;-' _,  veins
#        ;        `./   /`,    \-'
#        | right   /   |  ;\   |\
#        | atrium ;_,._|_,  `, ' \
#        |        \    \ `       `,
#        `      __ `    \   left  ;,
#         \   ,'  `      \,  ventricle
#          \_(            ;,      ;;
#          |  \           `;,     ;;
# inferior |  |`.          `;;,   ;'
#vena cava |  |  `-.        ;;;;,;'
#          |  |    |`-.._  ,;;;;;'
#          |  |    |   | ``';;;'  FL
#                  aorta

#Barve za cmap
N = 256
blueColorMatrix = np.ones((N, 4))
blueColorMatrix[:, 0] = np.linspace(40/N, 40/N, N)
blueColorMatrix[:, 1] = np.linspace(50/N, 50/N, N)
blueColorMatrix[:, 2] = np.linspace(114/N, 114/N, N)
blueColorMatrix[:, 3] = np.linspace(alpha_boost, 1, N)
blueCmap = ListedColormap(blueColorMatrix)

redColorMatrix = np.ones((N, 4))
redColorMatrix[:, 0] = np.linspace(127/N, 127/N, N)
redColorMatrix[:, 1] = np.linspace(23/N, 23/N, N)
redColorMatrix[:, 2] = np.linspace(16/N, 16/N, N)
redColorMatrix[:, 3] = np.linspace(alpha_boost, 1, N)
redCmap = ListedColormap(redColorMatrix)

if DataGroup2:
    DEBUG_txt.append(f' ( 3/3 ) Generating PROTOMAP results [{time.time()}]  ...')
    print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")
else:
    DEBUG_txt.append(f' ( 3/3 ) Generating SEQUENCE ALIGN results [{time.time()}] ...')
    print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")

if not os.path.isdir(work_dir):
    os.mkdir(work_dir)

if DataGroup2:
    sorting_column_caption = "K-factor"
else:
    sorting_column_caption = "coverage"

def drawPeptographs(protein):
    #SUMMARY
    total_coverage = [False] * len(FASTA[protein]) #ToDo: preverjanje če FASTA[protein] sploh obstaja
    #PRIPRAVIM FAJLE
    result_file = work_dir+slugify(protein.replace("|", "_"))+'.txt' #TXT
    result_figure = work_dir+slugify(protein.replace("|", "_"))+OUTPUT_EXTENTION #SVG / PNG
    if os.path.exists(result_file): #TXT
        os.remove(result_file)
    if os.path.exists(result_figure): #SVG
        os.remove(result_figure)
    
    if TXT_OUTPUT:
        result = codecs.open(result_file, "w", "utf-8") #TXT Pomoje to upočasnjuje celo skripto, zakaj ne bi tega zapisal na koncu??
        result.write(u'\ufeff')
        if protein in FASTA and len(FASTA[protein]) > 0:
            result.write("SEQ:\t"+FASTA[protein]+"\n") #ToDo: nastavit ustrezno število tabov
        else:
            print(f"{bcolors.FAIL}ERROR! Given fasta file doesn't match fasta results from MaxQuant!\nERROR found at entry \"{protein}\"!{bcolors.ENDC}")
            pause()
            exit()
        
    ###
    # IMAGE GENERATION START
    ###
    protein_length = len(FASTA[protein])
    labels = []
    if setMW:
        labels.append("") # TOP
        labels += [str(i) for i in setMW_Labels]
        labels.append("") # BOTTOM
    else:
        for i in range(N_max):
            labels.append(str(i+1))
    #graf
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [2, 1.2]})
    fig.suptitle(FASTA_name[protein])

    #leva stran
    ax1.set_xlim(1, protein_length)
    # UniProt Features
    if useUniProt_features:
        ax1.set_ylim(N_max+0.5+(uniProtFeature_coeficient*(N_max+1)), 0.5) # če features, potem k*(N_max+1) večji graf
        ax1.axhline(y=N_max+0.5, color='#000000', linestyle='-', linewidth=0.4)
        ###
        # Risanje featurjev
        ###
        legendbar_height =  uniProtFeature_coeficient*(N_max+1) / len(useUniProt_features)
        for i in range(len(useUniProt_features)):
            for feature in FASTA_features[protein][UniProt_Feature_types[useUniProt_features[i]]]:
                # feature = ['Signal peptide', 1, 16, False, '#00ff00']
                ax1.add_patch(Rectangle((feature[1], N_max+0.5+(i*legendbar_height)+(0.1*legendbar_height)), feature[2], 0.8*legendbar_height, alpha=UniProt_Features_alpha, linewidth=0.4, edgecolor='k', facecolor=str(feature[4])))
        
        # narišem še skupno intenziteto
        if DataGroup2: # če je samo ena skupina je res pointles in čudno zgleda
            total_height = uniProtFeature_coeficient*(N_max+1)
            ax2.barh(N_max+0.5+(total_height*0.3), Total_intensity_1[protein], height=(0.4*total_height), color='#283272', alpha=0.5)
            ax2.barh(N_max+0.5+(total_height*0.7), Total_intensity_2[protein], height=(0.4*total_height), color='#7f1710', alpha=0.5)
    else:
        ax1.set_ylim(N_max+0.5, 0.5)
       
    # y-os
    y_ticks = []
    if setMW:
        ax1.set_ylabel('Molecular weight [kDa]')
        for i in range(N_max+1):
            y_ticks.append(i+0.5)
    else:
        ax1.set_ylabel('Slice')
        for i in range(N_max):
            y_ticks.append(i+1)

    # obe x osi
    ax1.set_yticks(y_ticks, labels)
        
    ax1.set_xlabel('Amino acid position', labelpad=0)
    ax1.set_xticks([1, protein_length])
                
    secAx1 = ax1.twiny()
    secAx1.set_xlim(1, protein_length)
    secAx1.set_xticks([1, protein_length], ["N-term", "C-term"])
        
    #desna stran
    ax2.set_xlim(0, 1.01)
        
    if useUniProt_features:
        ax2.set_ylim(N_max+0.5+(uniProtFeature_coeficient*(N_max+1)), 0.5)
        ax2.axhline(y=N_max+0.5, color='#000000', linestyle='-', linewidth=0.4)
    else:
        ax2.set_ylim(N_max+0.5, 0.5)
        
    ax2.set_xlabel('Slice intensity', labelpad=10)
    ax2.tick_params(bottom=False)
    ax2.set_yticks(list(range(1, N_max+1)))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
        
    #horizontalne črte
    if horizontal_lines:
        for i in range(N_max):
            ax1.axhline(y=i+0.5, color='#000000', linestyle='-', linewidth=0.1)
            ax2.axhline(y=i+0.5, color='#000000', linestyle='-', linewidth=0.1)
        
    #nastavitev zamikov
    label_length_delta = len(max(labels, key=len)) - 5
    if label_length_delta > 0:
        y_offest = 0.1 + (label_length_delta * 0.009)
    else:
        y_offest = 0.1
    plt.subplots_adjust(left = y_offest, bottom = 0.11, right = 0.95, top = 0.87, wspace=0)
        
    #cmap
    norm = plt.Normalize(alpha_boost, 1)
    smBlue = plt.cm.ScalarMappable(cmap=blueCmap, norm=norm)
    smBlue.set_array([])
        
    if DataGroup2:
        smRed = plt.cm.ScalarMappable(cmap=redCmap, norm=norm)
        smRed.set_array([])
            
        fig.text(0.907-((y_offest-0.1)/5), 0.2, DataGroup1_name, rotation='vertical')
        fig.text(0.879-((y_offest-0.1)/5), 0.2, DataGroup2_name, rotation='vertical')
        BlueCbar = plt.colorbar(smBlue, ax=ax2, ticks=[alpha_boost, 1], pad=-0.05)
        #BlueCbar.set_ticklabels(['min', 'max'])
        BlueCbar.set_ticklabels([str(round(alpha_boost*100))+"%", str(100)+"%"])
        RedCbar = plt.colorbar(smRed, ax=ax2, ticks=[alpha_boost, 1], pad=0.05)
        RedCbar.set_ticks([])
    else:
        BlueCbar = plt.colorbar(smBlue, ax=ax2, ticks=[alpha_boost, 1])
        #BlueCbar.set_ticklabels(['min', 'max'])
        BlueCbar.set_ticklabels([str(round(alpha_boost*100))+"%", str(100)+"%"])
        
    #Začnem obdelavo posameznih slajsov
    slice = 1
    for i in range(N_max):
        #DEBUG_txt.append(f'NOTICE: working on slice {slice} in protein {protein} (DataGroup2={DataGroup2}) ...')
        if DataGroup2:
            experiment = data_titles[DataGroup1[i]]+'/'+data_titles[DataGroup2[i]]
        else:
            experiment = data_titles[DataGroup1[i]]
            
        coverage_list_Group1 = list(empty_char * len(FASTA[protein]))
        if DataGroup2:
            coverage_list_Group2 = list(empty_char * len(FASTA[protein]))
           
        for peptide in DB[protein]:
            peptide_start = PEPTIDES_qualitative[peptide][0]
            peptide_length = PEPTIDES_qualitative[peptide][1]
                
            try:
                peptide_intensity_Group1 = PEPTIDES_quantitative_avg_norm[peptide][i]
            except KeyError:
                DEBUG_txt.append(f'ERROR: Peptide "{peptide}" not found in variable PEPTIDES_quantitative_avg_norm[peptide][{i}]!')
                write_debug()
                print(bcolors.FAIL+DEBUG_txt[-1]+bcolors.ENDC)
                pause()
                peptide_intensity_Group1 = PEPTIDES_quantitative_avg_norm[peptide][i]
                   
             
            if DataGroup2:
                peptide_intensity_Group2 = PEPTIDES_quantitative_avg_norm[peptide][N_max+i]
                 
            alpha_intensity_Group1 = (1-alpha_boost)*peptide_intensity_Group1 + alpha_boost
            if DataGroup2:
                alpha_intensity_Group2 = (1-alpha_boost)*peptide_intensity_Group2 + alpha_boost
            #PROTOMAP
            if DataGroup2:
                #DataGroup1
                if peptide_intensity_Group1 > peptide_threshold:
                    # GRAPH(SVG) OUTPUT MAIN
                    ax1.add_patch(Rectangle((peptide_start, slice-(bar_height/2)), peptide_length, bar_height/2,alpha=alpha_intensity_Group1 , linewidth=0, color='#283272'))
                    # TXT RESULTS
                    if TXT_OUTPUT:
                        if peptide_start + peptide_length <= len(coverage_list_Group1):
                            for n in range(peptide_length):
                                total_coverage[peptide_start + n] = True
                                if n > 0 and n < (peptide_length-1):
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'mid')
                                elif n == 0:
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'start')
                                elif n == (peptide_length-1):
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'end')
                        else:
                            DEBUG_txt.append(f"WARNING: peptide \"{peptide}\" does not fit onto protein \"{protein}\"! Disregarding it in \"{experiment}\" ({PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1]} > {len(coverage_list_Group1)})")
                            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                elif peptide_intensity_Group1 > 0 and NOTICE:
                    DEBUG_txt.append(f"NOTICE: peptide \"{peptide}\" falls out of {round(peptide_threshold*100, 0)}% threshold ({round(peptide_intensity_Group1*100, 3)}%). Disregarding it in {protein} {data_titles[DataGroup1[i]]}")
                    print(DEBUG_txt[-1])
                #DataGroup2
                if peptide_intensity_Group2 > peptide_threshold:
                    # GRAPH(SVG) OUTPUT MAIN
                    ax1.add_patch(Rectangle((peptide_start, slice), peptide_length, bar_height/2,alpha=alpha_intensity_Group2 , linewidth=0, color='#7f1710'))
                    # TXT RESULTS
                    if TXT_OUTPUT:
                        if peptide_start + peptide_length <= len(coverage_list_Group2):
                            for n in range(peptide_length):
                                total_coverage[peptide_start + n] = True
                                if n > 0 and n < (peptide_length-1):
                                    coverage_list_Group2[peptide_start + n] = sequence_adjust(coverage_list_Group2[peptide_start + n], 'mid')
                                elif n == 0:
                                    coverage_list_Group2[peptide_start + n] = sequence_adjust(coverage_list_Group2[peptide_start + n], 'start')
                                elif n == (peptide_length-1):
                                    coverage_list_Group2[peptide_start + n] = sequence_adjust(coverage_list_Group2[peptide_start + n], 'end')
                        else:
                            DEBUG_txt.append(f"WARNING! Peptide \"{peptide}\" does not fit onto protein \"{protein}\"! Disregarding it in \"{experiment}\" ({PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1]} > {len(coverage_list_Group2)})")
                            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                elif peptide_intensity_Group2 > 0 and NOTICE:
                    DEBUG_txt.append(f"NOTICE: peptide \"{peptide}\" falls out of {round(peptide_threshold*100, 0)}% threshold ({round(peptide_intensity_Group2*100, 3)}%). Disregarding it in {protein} {data_titles[DataGroup2[i]]}")
                    print(DEBUG_txt[-1])
            #PEPTIDEALIGN
            else:
                #DataGroup1
                if peptide_intensity_Group1 > peptide_threshold:
                    # GRAPH(SVG) OUTPUT MAIN
                    ax1.add_patch(Rectangle((peptide_start, slice-(bar_height/2)), peptide_length, bar_height, alpha=alpha_intensity_Group1 , linewidth=0, color='#283272'))
                    # TXT RESULTS
                    if TXT_OUTPUT:
                        if peptide_start + peptide_length <= len(coverage_list_Group1):
                            for n in range(peptide_length):
                                total_coverage[peptide_start + n] = True
                                if n > 0 and n < (peptide_length-1):
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'mid')
                                elif n == 0:
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'start')
                                elif n == (peptide_length-1):
                                    coverage_list_Group1[peptide_start + n] = sequence_adjust(coverage_list_Group1[peptide_start + n], 'end')
                        else:
                            DEBUG_txt.append(f"WARNING! Peptide \"{peptide}\" does not fit onto protein \"{protein}\"! Disregarding it in \"{experiment}\" ({PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1]} > {len(coverage_list_Group1)})")
                            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                elif peptide_intensity_Group1 > 0 and NOTICE:
                    DEBUG_txt.append(f"NOTICE: peptide \"{peptide}\" falls out of {round(peptide_threshold*100, 0)}% threshold ({round(peptide_intensity_Group1*100, 3)}%). Disregarding it in {protein} {data_titles[DataGroup1[i]]}")
                    print(DEBUG_txt[-1])
        # dodam posamezni slice v txt in narišem skupno intenziteto slica + napako, če je podana
        if TXT_OUTPUT:
            coverage_Group1 = ''.join(coverage_list_Group1)
            result.write(data_titles[DataGroup1[i]]+"\t"+coverage_Group1+"\n")
            del coverage_list_Group1
        if DataGroup2:
            if TXT_OUTPUT:
                coverage_Group2 = ''.join(coverage_list_Group2)
                result.write(data_titles[DataGroup2[i]]+"\t"+coverage_Group2+"\n")
                del coverage_list_Group2
            # rišem dve intenziteti
            if error_bar:
                this_err_1 = SLICE_quantitative_err[protein][slice]
                this_err_2 = SLICE_quantitative_err[protein][N_max+slice]
                err_height = 2
            else:
                err_height = this_err_1 = this_err_2 = 0
              
            ax2.barh(slice-(bar_width/4), SLICE_quantitative[protein][slice], xerr=this_err_1, height=(bar_width/2), color='#283272', alpha=0.5, ecolor='black', capsize=err_height)
            ax2.barh(slice+(bar_width/4), SLICE_quantitative[protein][N_max+slice], xerr=this_err_2, height=(bar_width/2), color='#7f1710', alpha=0.5, ecolor='black', capsize=err_height)

        else:
            #rišem eno intenziteto
            if error_bar:
                this_err = SLICE_quantitative_err[protein][slice]
                err_height = 2
            else:
                err_height = this_err = 0
            ax2.barh(slice, SLICE_quantitative[protein][slice], xerr=this_err, height=bar_width, color='#283272', alpha=0.5, ecolor='black', capsize=err_height)
         
        slice+=1
      
    #SUMMARY
    SUMMARY[protein].append(round((sum(total_coverage) / len(total_coverage))*100, 1)) #total coverage
    if DataGroup2:
        SUMMARY[protein].append(round(K_factor[protein]*(sum(total_coverage) / len(total_coverage)), 2)) # K-factor - morda ga bi lahko množil s pokritostjo, da ne bo polno proteinov, ki so samo malo zaznani?
    else:
        SUMMARY[protein].append(SUMMARY[protein][-1]) # če ni dataset2, potem izpisuje coverage
    
    #TXT KONEC
    if TXT_OUTPUT:
        result.close()
    #SVG KONEC
    #plt.savefig(result_figure, format='svg', dpi=1200)
    if SVG_OUTPUT:
        plt.savefig(result_figure, format=OUTPUT_EXTENTION[1:], dpi=1200, transparent=transparent_bg)
    else:
        plt.savefig(result_figure, format=OUTPUT_EXTENTION[1:], dpi=150, transparent=transparent_bg)
    plt.close()
    return protein

# MULTIPROCESS
if Multiproc and N_write_threads > 1:
    with alive_bar(len(multi_IDs), dual_line=True, theme='classic') as bar:
        results_draw = ThreadPool(N_write_threads).imap_unordered(drawPeptographs, multi_IDs)
        for r in results_draw:
            bar.text = f'Drawn peptograf for {r}'
            bar()
else:
    with alive_bar(len(multi_IDs), dual_line=True, theme='classic') as bar:
        for id in multi_IDs:
            drawPeptographs(id)
            bar.text = f'Drawn peptograf for {id}'
            bar()
   
###
# Risanje legend.svg
###
if ONLINE and useUniProt_features:
    figl, axl = plt.subplots(figsize=(6,4))
    axl.axis(False)

    for feature_ID in useUniProt_features:
        if feature_ID == 0: #molecular processing
            THIS_LEGEND = LEGEND_Processing
            lloc = "upper center"
        elif feature_ID == 1:
            THIS_LEGEND = LEGEND_Region
            lloc = "center"
        elif feature_ID == 2:
            THIS_LEGEND = LEGEND_2ndary
            lloc = "lower center"
        else:
            THIS_LEGEND = False
            DEBUG_txt.append(f"ERROR! Unknown feature ID {feature_ID}! Skipping legend creation for unknown feature...")
        
        # Generiranje THIS_LEGEND
        if THIS_LEGEND:
            #UniProt_Feature_types[feature_ID]
            legend_labels = legend_handles = []
            for i in range(len(THIS_LEGEND)):
                legend_labels.append(THIS_LEGEND[i][0])
            legend_handles = [plt.Rectangle((0,0),1,1, color=THIS_LEGEND[i][1], alpha=UniProt_Features_alpha) for i in range(len(THIS_LEGEND))]
            axl = plt.gca().add_artist(plt.legend(legend_handles, legend_labels, loc=lloc, ncol=3, title=UniProt_Feature_types[feature_ID]))
        
        del THIS_LEGEND
    figl.savefig(f'{work_dir}legend{OUTPUT_EXTENTION}')

###
##### Generiranje HTML datoteke z rezultati
###
DEBUG_txt.append(f'ALL DONE! Generating {exp_name}.html results file in "{os.path.dirname(peptides_file)}" ... ')
print("\n"+DEBUG_txt[-1])
filenames = []
#for protein,summary_data in SUMMARY.items():
for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][1], reverse=True):
    filenames.append(summary_data[0])

groups_summary = "DataGroup1: "
for key in DataGroup1:
    groups_summary += data_titles[key]+", "

groups_summary = groups_summary[0:-2]+"<br>"
if(DataGroup2):
    groups_summary += "DataGroup2: "
    for key in DataGroup2:
        groups_summary += data_titles[key]+", "
    groups_summary = groups_summary[0:-2]+"<br>Group normalization?"+str(group_normalization)

HTML_CODE = '''
<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
	<head>
		<meta http-equiv="content-type" content="text/html; charset=UTF-8" />
		<style type="text/css">
			body {margin:0;padding:0;height:100%;}
			div#list {clear:both;position:fixed;top:0;height:25%;width:100%;background-color:#9999cc;color:#FFFFFF;font-size:16px;overflow-y: scroll;}
			div#content {clear:both;position:fixed;bottom:0;height:75%;width:100%;background-color:#dddddd;color:#000000;font-size:16px;overflow-y:scroll;text-align:center;}
			.main_colum {width:50%;padding-left:2em;}
			.sub_colum {width:12.5%;padding-left:2em;}
			th {cursor: pointer;}
		</style>
		<script>
			const PROTEINS = '''+json.dumps(filenames)+''';
			N = PROTEINS.length - 1;
			last_call = -1;
			
			function loadProtein(ID) {
				document.getElementById("content").innerHTML = '<img src="'''+exp_name+"""/'+PROTEINS[ID]+'"""+OUTPUT_EXTENTION+'''" alt="'+PROTEINS[ID]+'"><br><img src="'''+exp_name+'''/legend'''+OUTPUT_EXTENTION+'''">';
				last_call = ID;
			}
			
			document.addEventListener('keydown', function(e) {
				switch (e.keyCode) {
					case 37:
					case 38:
						if(last_call > 0) {
							loadProtein(last_call-1)
						}/* else {
							loadProtein(N)
						}*/
						break;
					case 39:
					case 40:
						if(last_call < N || last_call == -1) {
							loadProtein(last_call+1)
						}/* else {
							loadProtein(0)
						}*/
						break;
				}
			});
		</script>
	</head>
	<body>
		<div id="list">
			<table style="width:100%">
				<thead>
					<tr style="text-align:left; background-color:#aaaadd;">
						<th class="main_colum">Protein</th>
                        <th class="sub_colum">'''+sorting_column_caption+'''</th>
                        <th class="sub_colum">% coverage</th>
						<th class="sub_colum"># peptides</th>
						<th class="sub_colum">Gross intensity</th>
					</tr>
				</thead>
				<tbody>
'''
n = 0
#for protein,summary_data in SUMMARY.items():
if useUniProt:
    for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][1], reverse=True):
        HTML_CODE += '''				    <tr onclick="loadProtein('''+str(n)+''')">
						<td>'''+FASTA_name[protein]+' (UniProt ID: '+protein+')'+'''</td>
                        <td>'''+str(summary_data[4])+'''</td>
                        <td>'''+str(summary_data[3])+'''</td>
						<td>'''+str(summary_data[2])+'''</td>
						<td>'''+str(summary_data[1])+'''</td>
                  </tr>
'''
        n+=1
else:
    for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][1], reverse=True):
        HTML_CODE += '''				    <tr onclick="loadProtein('''+str(n)+''')">
						<td>'''+protein+'''</td>
                        <td>'''+str(summary_data[4])+'''</td>
                        <td>'''+str(summary_data[3])+'''</td>
						<td>'''+str(summary_data[2])+'''</td>
						<td>'''+str(summary_data[1])+'''</td>
                  </tr>
'''
        n+=1

HTML_CODE += '''
				</tbody>
			</table>
            <script>
                // Sorting .. addapted from Nick Grealy
				const getCellValue = (tr, idx) => tr.children[idx].innerText || tr.children[idx].textContent;

                const comparer = (idx, asc) => (a, b) => ((v1, v2) => 
                    v1 !== '' && v2 !== '' && !isNaN(v1) && !isNaN(v2) ? v1 - v2 : v1.toString().localeCompare(v2)
                    )(getCellValue(asc ? a : b, idx), getCellValue(asc ? b : a, idx));

                document.querySelectorAll('th').forEach(th => th.addEventListener('click', (() => {
                    const table = th.closest('table');
                    const tbody = table.querySelector('tbody');
                        Array.from(tbody.querySelectorAll('tr'))
                            .sort(comparer(Array.from(th.parentNode.children).indexOf(th), this.asc = !this.asc))
                            .forEach(tr => tbody.appendChild(tr) );
			})));
			</script>
            </script>
		</div>
		<div id="content"><br>'''+"Experiment \""+usr_name+"\"<br>LFQ intensities used? "+str(use_LFQ)+"<br>Relative threshold: "+str((peptide_threshold * 100))+"%<br>"+groups_summary
if useUniProt:
    HTML_CODE += '<br><img src="'+exp_name+'/legend'+OUTPUT_EXTENTION+'">'
HTML_CODE += '<br><br>DEBUG TXT<br>'+'<br>'.join(DEBUG_txt)+'''</div>
	</body>
</html>
'''

HTML = codecs.open(os.path.dirname(peptides_file)+'/'+exp_name+".html", "w", "utf-8") #TXT
HTML.write(u'\ufeff')
HTML.write(HTML_CODE)
HTML.close()

###
##### Konec + DEBUG
###
process_time = time.time() - continue_time
if DEBUG:
    print("")
    write_debug()

DEBUG_txt.append(f"DONE! [{time.time()}]")
print(f"\n{bcolors.OKGREEN}{DEBUG_txt[-1]}{bcolors.ENDC}")
DEBUG_txt.append(f"Read {line_count} peptides in {round(read_time)*1000} ms and processed results in {round(process_time/60)} minutes.")
print(f"{bcolors.OKBLUE}{DEBUG_txt[-1]}{bcolors.ENDC}\n")
pause()