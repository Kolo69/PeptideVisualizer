#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# (c) 2022
#     Matej Kolarič (matej.kolaric@ijs.si), Robert Vidmar, Marko Fonović
#     Biochemistry, molecular and structural biology - B1, Jozef Stefan Institute, Jamova cesta 39, Ljubljana Slovenia
#     under GPLv3 (http://www.gnu.org/licenses/gpl-3.0.html)
#
# Contributors: Tilen Sever, Marija Grozdanić, Sara Ivanovski, Andreja Kozak
#
#

##############################################################################
#               Nastavitve
##############################################################################

#ToDo: implement no CLI option if GUI available and rather provide CLI after crash
DEBUG = False
NOTICE = False
SELFUPDATE = True
peptograph_exp_spacing = 2 # v odstotkih
bar_height = 0.9
bar_width = 0.8
alpha_boost = 0.05
default_peptide_threshold = 5 # v odstotkih
uniProtFeature_coeficient = 0.07 # +x% višine grafa je legenda
UniProt_Features_alpha = 0.5
err_width = 3
err_cap = 3
err_capthick = 2
replicates_MAX = 10000
transparent_bg = False
TXT_OUTPUT = True
AA_PER_LINE = 100
AA_STEP = 10
Multiproc = True
N_download_threads = 16 # včasih zlaufa lepo z 100+, ampak obstaja verjetnost da se sesuje
N_write_threads = 1 # ČE >1 VČASI POSTANE RNG!

GUI = True
DARK_MODE = True
GUI_list_width = 26
GUI_list_height = 21
GUI_label_width = GUI_list_width-6
GUI_checkbox_width = GUI_list_width-9
GUI_list_pad = 5
GUI_std_pad = 10
GUI_checkbox_lvl2_pad = 15
GUI_number_pad = 90
GUI_button_pad = 50
DEV_DPI = 72.0
DEV_scale = 1.0

PEPTOGRAPH_COLORS_PRESETS = {
    "Diadic": ["283272", "7f1710"],
    "Triadic": ["00a88f", "a64499", "f78f1e"],
    "Tetradic": ["0a3aab", "ce0073", "9aee00", "ffaa00"],
    "Hexadic": ["cc0098", "670099", "ffcc00", "ff6600", "66cc00", "0ab4c3"]
}

UniProt_Feature_types = ["Processing", "Region", "Secondary structure"]
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
VERSION = '1.7'

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

Contributors: Tilen Sever, Marija Grozdanić, Sara Ivanovski, Andreja Kozak

##########################################################################################################"""

#globalne spremenljivke
usr_name = ""
DataGroup1_name = ""
DataGroup2_name = ""
DataGroup_names = []
DataGroup_values = []
EXP_N = 0
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
MISSING_DEPENDENCIES = False

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

if OS == "MacOS" or OS == "Linux":
    if os.geteuid() != 0:
        missing_root_privileges = True
        DEBUG_txt.append(f"WARNING: No root privileges! The script will fail if any dependencies are missing! os.geteuid() = {os.geteuid()}")
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
    if OS == "Linux":
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
    MISSING_DEPENDENCIES = True
    DEBUG_txt.append(f'WARNING! Python module numpy not found!')
    if missing_root_privileges:
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
        pause()
        exit()
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install numpy\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling numpy ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'numpy', "--quiet"])
    time.sleep(1)
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed numpy!")
    print(DEBUG_txt[-1]+"\n")

# REQUIRED: matplotlib
DEBUG_txt.append('NOTICE: Importing matplotlib ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    import matplotlib
except ImportError:
    MISSING_DEPENDENCIES = True
    DEBUG_txt.append(f'WARNING! Python module matplotlib not found!')
    if missing_root_privileges:
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
        pause()
        exit()
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install matplotlib\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling matplotlib ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'matplotlib', "--quiet"])
    time.sleep(1)
    #import matplotlib
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed matplotlib!")
    print(DEBUG_txt[-1]+"\n")

if not MISSING_DEPENDENCIES:
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
    MISSING_DEPENDENCIES = True
    DEBUG_txt.append(f'WARNING! Python module requests not found!')
    if missing_root_privileges:
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
        pause()
        exit()
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install requests\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling requests ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'requests', "--quiet"])
    time.sleep(1)
    #import requests
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed requests!")
    print(DEBUG_txt[-1]+"\n")

# datumi za apdejt
DEBUG_txt.append('NOTICE: Importing dateutil.parser ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    from dateutil.parser import parse as parsedate
except ImportError:
    MISSING_DEPENDENCIES = True
    DEBUG_txt.append(f"WARNING! Python module dateutil not found!")
    if missing_root_privileges:
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
        pause()
        exit()
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install python-dateutil\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling python-dateutil ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'python-dateutil', "--quiet"])
    time.sleep(1)
    #from dateutil.parser import parse as parsedate
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed dateutil!")
    print(DEBUG_txt[-1]+"\n")

# REQUIRED: alive_progress
DEBUG_txt.append('NOTICE: Importing alive_progress ...')
if NOTICE:
    print(DEBUG_txt[-1])
try:
    from alive_progress import alive_bar
except ImportError:
    MISSING_DEPENDENCIES = True
    DEBUG_txt.append(f"WARNING! Python module alive-progress not found!")
    if missing_root_privileges:
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
        pause()
        exit()
    print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install alive-progress\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
    pause()
    print("Istalling alive-progress ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'alive-progress', "--quiet"])
    time.sleep(1)
    #from alive_progress import alive_bar
    DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed alive_progress! Please restart the script.")
    print(DEBUG_txt[-1]+"\n")

###
# GUI (tkinter) import + basic checks
###
if GUI:
    DEBUG_txt.append('NOTICE: Importing tkinter ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        import tkinter
    except ImportError:
        GUI = False
        DEBUG_txt.append(f"WARNING! Python module tkinter (python3-tk) not found!")
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} is unable to install it for you.{bcolors.ENDC}")
        if OS == "Linux":
            print(f"Please run \"sudo apt-get install python3-tk\" or equivalent to use python3-tk ...")
        elif OS == "MacOS":
            print(f"Please run \"brew install python-tk\" to use python3-tk ...")
        else:
            print(f"I have no idea how tkinter is not installed on OS={OS} ... Please try updating python to use GUI.")
    
    DEBUG_txt.append('NOTICE: Importing askopenfilename ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        from tkinter.filedialog import askopenfilename
        from tkinter.colorchooser import askcolor
    except:
        GUI = False
    
    DEBUG_txt.append('NOTICE: Importing FigureCanvasTkAgg ...')
    if NOTICE:
        print(DEBUG_txt[-1])
    try:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #for MW setting
    except:
        MISSING_DEPENDENCIES = True
        DEBUG_txt.append(f"WARNING! Python module python3-pil.imagetk not found!")
        if missing_root_privileges:
            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}\n{bcolors.ERROR}ERROR! insufficient privileges! Please run {PROGRAM_NAME} as sudoer or root! Terminating... {bcolors.ENDC}")
            pause()
            exit()
        print(f"{bcolors.WARNING}{DEBUG_txt[-1]} {PROGRAM_NAME} will try to install it for you. If the script fails run \"pip install pillow\" in command promt.{bcolors.ENDC}\nNOTICE: {PROGRAM_NAME} will fail without internet connection!")
        pause()
        print("Istalling python3-pil.imagetk ... ")
        subprocess.check_call([sys.executable, "-m", "pip", "install", 'pillow', "--quiet"])
        time.sleep(1)
        #from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        DEBUG_txt.append(f"{PROGRAM_NAME} successfully installed python3-pil.imagetk!")
        print(DEBUG_txt[-1]+"\n")
if GUI:
    try:
        import tkinter as tk
    except ImportError: # je to sploh modžno?
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
    try:
        tktest = tk.Tk()
        screen_w = tktest.winfo_screenwidth()
        screen_h = tktest.winfo_screenheight()
        screen_dpi = tktest.winfo_fpixels('1i')
        OS_scale = screen_dpi / 96 # 96 is a scale factor for 100%
        scale_f = screen_dpi / DEV_DPI
        screen_w_eff = round(screen_w / OS_scale)
        screen_h_eff = round(screen_h / OS_scale)
        wh = round(900 * OS_scale)
        ww = round(550 * OS_scale)
        tktest.withdraw()
        tktest.update()
        tktest.destroy()
        DEBUG_txt.append(f"NOTICE: resolution={screen_w}x{screen_h} dpi={screen_dpi} OS_scale={OS_scale} ({screen_w_eff}x{screen_h_eff})")
        DEBUG_txt.append(f"NOTICE: scale_f={scale_f} => {wh}x{ww}")
        if NOTICE:
            print(DEBUG_txt[-2]+"\n"+DEBUG_txt[-1])
    except:
        GUI = False

DEBUG_txt.append(f'NOTICE: Graphical User Interface = {GUI}')
if NOTICE:
    print(DEBUG_txt[-1])

if MISSING_DEPENDENCIES:
    DEBUG_txt.append(f'{bcolors.WARNING}WARNING: MISSING_DEPENDENCIES = {MISSING_DEPENDENCIES}! {PROGRAM_NAME} installed missing dependencies. Please restart the script...{bcolors.ENDC}')
    print(DEBUG_txt[-1])
    pause()
    exit()

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

# COLORS
PEPTOGRAPH_COLORS = PEPTOGRAPH_COLORS_PRESETS[list(PEPTOGRAPH_COLORS_PRESETS.keys())[0]]*100
class PeptideVisualizerGUI(tk.Tk):
    def __init__(self, *args, **kwargs):
        global ONLINE
        tk.Tk.__init__(self, *args, **kwargs)
        
        # DARK/LIGHT Mode
        if DARK_MODE:
            self.bg, self.fg, self.ac1, self.ac2 = ('#282828', 'white', '#404040', '#B3B3B3')
        else:
            self.bg, self.fg, self.ac1, self.ac2 = ('#f0f0ed', 'black', '#ffffff', '#E9DAC1')
        
        #self.tk.call('tk', 'scaling', scale_f)
        #self.font = ('Serif', round(12*scale_f))
        self.configure(bg=self.bg, borderwidth=0)
        self.title(PROGRAM_NAME)
        #self.geometry(f"{wh}x{ww}")
        self.resizable() #self.resizable(width=False, height=False)
        
        # Variables
        self.Experiment_name = tk.StringVar()
        self.replicates_var = tk.IntVar()
        self.replicates_var.set(1)
        self.peptide_threshold_var = tk.DoubleVar()
        self.peptide_threshold_var.set(default_peptide_threshold)
        self.saveSVG = tk.BooleanVar(value=True)
        
        self.DataSets_vars = tk.StringVar()
        
         # Imena eksperimentov
        self.DataGroup1_name = tk.StringVar()
        self.DataGroup2_name = tk.StringVar()
        self.DataGroup3_name = tk.StringVar()
        self.DataGroup4_name = tk.StringVar()
        self.DataGroup5_name = tk.StringVar()
        self.DataGroup6_name = tk.StringVar()
        self.DataGroup7_name = tk.StringVar()
        self.DataGroup8_name = tk.StringVar()
        self.DataGroup1_name.set("Dataset 1")
        self.DataGroup2_name.set("Dataset 2")
        self.DataGroup3_name.set("Dataset 3")
        self.DataGroup4_name.set("Dataset 4")
        self.DataGroup5_name.set("Dataset 5")
        self.DataGroup6_name.set("Dataset 6")
        self.DataGroup7_name.set("Dataset 7")
        self.DataGroup8_name.set("Dataset 8")
        # Spremenljivke s seznami
        self.DataGroup1_vars = tk.StringVar()
        self.DataGroup2_vars = tk.StringVar()
        self.DataGroup3_vars = tk.StringVar()
        self.DataGroup4_vars = tk.StringVar()
        self.DataGroup5_vars = tk.StringVar()
        self.DataGroup6_vars = tk.StringVar()
        self.DataGroup7_vars = tk.StringVar()
        self.DataGroup8_vars = tk.StringVar()
        
        # Nastavitve PeptideVisualizerja
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
        self.Experiment_label = tk.Label(text="Experiment Name:", width=GUI_label_width, anchor="e")
        self.Experiment_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.Experiment_label.grid(row=0, column=0, sticky="e")
        
        self.Experiment_text = tk.Entry(textvariable = self.Experiment_name,width=GUI_list_width)
        self.Experiment_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.Experiment_text.grid(row=0, column=1, pady=(10, 10))
        
        self.DataSets_name = tk.Label(text="Data from peptides.txt", width=GUI_label_width, anchor="w")
        self.DataSets_name.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.DataSets_name.grid(row=1, column=0, sticky="w", padx=(GUI_std_pad, 0))
        
        # Color presets
        color_options = []
        for preset_name in PEPTOGRAPH_COLORS_PRESETS.keys():
            color_options.append(preset_name)
        if len(color_options) != 0:
            self.preset_name = tk.StringVar(value="Color presets")
            self.preset_menu = tk.OptionMenu(self, self.preset_name, *color_options, command=self.set_color_preset)
            self.preset_menu.config(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, borderwidth=0, width=GUI_checkbox_width)
            self.preset_menu["menu"].config(bg=self.bg, fg=self.fg)
            self.preset_menu.grid(row=0, column=2, pady=(10, 10), padx=(GUI_label_width-GUI_checkbox_width+GUI_std_pad,GUI_std_pad))
        
        # Imena + barve
        self.DataGroup1_text = tk.Entry(textvariable = self.DataGroup1_name, width=GUI_label_width)
        self.DataGroup1_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup1_text.grid(row=1, column=1, sticky='W', padx=(GUI_list_pad,0)) # Samo en je po defaultu viden
        self.DataGroup1_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[0], command=self.chg1, font=('Helvetica', 6))
        self.DataGroup1_color.grid(row=1, column=1, sticky='E')
        
        self.DataGroup2_text = tk.Entry(textvariable = self.DataGroup2_name, width=GUI_label_width)
        self.DataGroup2_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup2_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[1], command=self.chg2, font=('Helvetica', 6))
        self.DataGroup3_text = tk.Entry(textvariable = self.DataGroup3_name,width=GUI_label_width)
        self.DataGroup3_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup3_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[2], command=self.chg3, font=('Helvetica', 6))
        self.DataGroup4_text = tk.Entry(textvariable = self.DataGroup4_name, width=GUI_label_width)
        self.DataGroup4_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup4_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[3], command=self.chg4, font=('Helvetica', 6))
        self.DataGroup5_text = tk.Entry(textvariable = self.DataGroup5_name, width=GUI_label_width)
        self.DataGroup5_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup5_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[4], command=self.chg5, font=('Helvetica', 6))
        self.DataGroup6_text = tk.Entry(textvariable = self.DataGroup6_name, width=GUI_label_width)
        self.DataGroup6_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup6_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[5], command=self.chg6, font=('Helvetica', 6))
        self.DataGroup7_text = tk.Entry(textvariable = self.DataGroup7_name, width=GUI_label_width)
        self.DataGroup7_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup7_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[6], command=self.chg7, font=('Helvetica', 6))
        self.DataGroup8_text = tk.Entry(textvariable = self.DataGroup8_name, width=GUI_label_width)
        self.DataGroup8_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup8_color = tk.Button(text='    ', bg='#'+PEPTOGRAPH_COLORS[7], command=self.chg8, font=('Helvetica', 6))
        
        # DataSet retrived from peptides.txt
        self.DataSets_list = tk.Listbox(listvariable=self.DataSets_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataSets_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataSets_list.grid(row=2, column=0, padx=(GUI_list_pad,0))
        self.DataSets_vars.set(value=data_titles_INT)
        
        # DataGroup1 (onlyone shown by default)
        self.DataGroup1_list = tk.Listbox(listvariable=self.DataGroup1_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup1_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.DataGroup1_list.grid(row=2, column=1, padx=(GUI_list_pad,0))
        self.add1 = tk.Button(text=' ADD ', command=self.add_1)
        self.add1.configure(bg=self.ac1, fg=self.fg)
        self.add1.grid(row=3, column=1, padx=(0,GUI_button_pad))
        self.remove1 = tk.Button(text=' DEL ', command=self.del_1)
        self.remove1.configure(bg=self.ac1, fg=self.fg)
        self.remove1.grid(row=3, column=1, padx=(GUI_button_pad,0))
        # DataGroup2
        self.DataGroup2_list = tk.Listbox(listvariable=self.DataGroup2_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup2_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add2 = tk.Button(text=' ADD ', command=self.add_2)
        self.add2.configure(bg=self.ac1, fg=self.fg)
        self.remove2 = tk.Button(text=' DEL ', command=self.del_2)
        self.remove2.configure(bg=self.ac1, fg=self.fg)
        # DataGroup3
        self.DataGroup3_list = tk.Listbox(listvariable=self.DataGroup3_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup3_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add3 = tk.Button(text=' ADD ', command=self.add_3)
        self.add3.configure(bg=self.ac1, fg=self.fg)
        self.remove3 = tk.Button(text=' DEL ', command=self.del_3)
        self.remove3.configure(bg=self.ac1, fg=self.fg)
        # DataGroup4
        self.DataGroup4_list = tk.Listbox(listvariable=self.DataGroup4_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup4_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add4 = tk.Button(text=' ADD ', command=self.add_4)
        self.add4.configure(bg=self.ac1, fg=self.fg)
        self.remove4 = tk.Button(text=' DEL ', command=self.del_4)
        self.remove4.configure(bg=self.ac1, fg=self.fg)
        # DataGroup5
        self.DataGroup5_list = tk.Listbox(listvariable=self.DataGroup5_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup5_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add5 = tk.Button(text=' ADD ', command=self.add_5)
        self.add5.configure(bg=self.ac1, fg=self.fg)
        self.remove5 = tk.Button(text=' DEL ', command=self.del_5)
        self.remove5.configure(bg=self.ac1, fg=self.fg)
        # DataGroup6
        self.DataGroup6_list = tk.Listbox(listvariable=self.DataGroup6_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup6_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add6 = tk.Button(text=' ADD ', command=self.add_6)
        self.add6.configure(bg=self.ac1, fg=self.fg)
        self.remove6 = tk.Button(text=' DEL ', command=self.del_6)
        self.remove6.configure(bg=self.ac1, fg=self.fg)
        # DataGroup7
        self.DataGroup7_list = tk.Listbox(listvariable=self.DataGroup7_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup7_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add7 = tk.Button(text=' ADD ', command=self.add_7)
        self.add7.configure(bg=self.ac1, fg=self.fg)
        self.remove7 = tk.Button(text=' DEL ', command=self.del_7)
        self.remove7.configure(bg=self.ac1, fg=self.fg)
        # DataGroup8
        self.DataGroup8_list = tk.Listbox(listvariable=self.DataGroup8_vars, height=GUI_list_height, width=GUI_list_width, selectmode='extended')
        self.DataGroup8_list.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.add8 = tk.Button(text=' ADD ', command=self.add_8)
        self.add8.configure(bg=self.ac1, fg=self.fg)
        self.remove8 = tk.Button(text=' DEL ', command=self.del_8)
        self.remove8.configure(bg=self.ac1, fg=self.fg)
        
        # Gumb za dodat nove eksperimente
        self.AddColumn = tk.Button(text="  + EXP  ", command=self.another_exp_pls)
        self.AddColumn.configure(bg=self.ac1, fg=self.fg)
        self.AddColumn.grid(row=2, column=2, padx=(GUI_list_pad,0))
        
        # Everything else - first column
        self.useLFQ_check = tk.Checkbutton(text="use LFQ", variable=self.useLFQ, width=GUI_checkbox_width, onvalue=True, offvalue=False, anchor="w")
        self.useLFQ_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useLFQ_check.grid(row=4, column=0, sticky='w', padx=(GUI_std_pad,0))
        
        self.replicates_label = tk.Label(text="Replicates  N = ", width=GUI_label_width, anchor="w")
        self.replicates_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.replicates_label.grid(row=5, column=0, sticky='w', padx=(GUI_std_pad,0))
        self.replicates_var_text = tk.Entry(textvariable = self.replicates_var, width=3)
        self.replicates_var_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.replicates_var_text.grid(row=5, column=0, sticky='w', padx=((GUI_number_pad+6)*OS_scale,0))# +6 ker je text daljši en znak
        
        self.peptide_threshold_label = tk.Label(text="Threshold  T =         %", width=GUI_label_width, anchor="w")
        self.peptide_threshold_label.configure(bg=self.bg, fg=self.fg, activebackground=self.bg, borderwidth=0)
        self.peptide_threshold_label.grid(row=6, column=0, sticky='w', padx=(GUI_std_pad,0))
        self.peptide_threshold_var_text = tk.Entry(textvariable = self.peptide_threshold_var, width=3)
        self.peptide_threshold_var_text.configure(bg=self.ac1, fg=self.fg, borderwidth=0)
        self.peptide_threshold_var_text.grid(row=6, column=0, sticky='w', padx=(GUI_number_pad*OS_scale,0))
        
        self.saveSVG_check = tk.Checkbutton(text="Save as .svg", variable=self.saveSVG, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)
        self.saveSVG_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.saveSVG_check.grid(row=7, column=0, sticky='w', padx=(GUI_std_pad,0))
        
        # second column
        self.horizontal_lines_check = tk.Checkbutton(text="Lines between slices", variable=self.horizontal_lines, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)
        self.horizontal_lines_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.horizontal_lines_check.grid(row=4, column=1, sticky='w', padx=(GUI_std_pad,0))
        
        self.setMW_check = tk.Checkbutton(text="Molar mass on Y-ax", variable=self.setMW, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)
        self.setMW_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.setMW_check.grid(row=5, column=1, sticky='w', padx=(GUI_std_pad,0))
        
        self.error_bar_check = tk.Checkbutton(text="Display error bars", variable=self.error_bar, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)#, state="disabled")
        self.error_bar_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.error_bar_check.grid(row=6, column=1, sticky='w', padx=(GUI_std_pad,0))
        
        self.group_normalization_check = tk.Checkbutton(text="Group normalization", variable=self.group_normalization, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False, command=self.group_normalization_alert)
        self.group_normalization_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.group_normalization_check.grid(row=7, column=1, sticky='w', padx=(GUI_std_pad,0))
        
        #third column
        self.useUniProt_check = tk.Checkbutton(text="use UniProt online DB", variable=self.useUniProt, width=GUI_checkbox_width, onvalue=True, offvalue=False, command=self.chgUniProt_check)
        self.useUniProt_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures1_check = tk.Checkbutton(text=UniProt_Feature_types[0], variable=self.useFeatures1, width=GUI_checkbox_width, anchor='w', onvalue=True, offvalue=False)
        self.useFeatures1_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures2_check = tk.Checkbutton(text=UniProt_Feature_types[1], variable=self.useFeatures2, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)
        self.useFeatures2_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        self.useFeatures3_check = tk.Checkbutton(text=UniProt_Feature_types[2], variable=self.useFeatures3, width=GUI_checkbox_width, anchor="w", onvalue=True, offvalue=False)
        self.useFeatures3_check.configure(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, selectcolor=self.ac1, borderwidth=0)
        
        if not ONLINE:
            self.useUniProt_check.config(text="UniProt unavailable (you are offline)", state="disabled")
            self.useFeatures1_check.config(state="disabled")
            self.useFeatures2_check.config(state="disabled")
            self.useFeatures3_check.config(state="disabled")  
        self.useUniProt_check.grid(row=4, column=2, sticky='w', padx=(GUI_std_pad,0))
        self.useFeatures1_check.grid(row=5, column=2, sticky='w', padx=(GUI_checkbox_lvl2_pad,0))
        self.useFeatures2_check.grid(row=6, column=2, sticky='w', padx=(GUI_checkbox_lvl2_pad,0))
        self.useFeatures3_check.grid(row=7, column=2, sticky='w', padx=(GUI_checkbox_lvl2_pad,0))
        
        self.VisualizerDo_btn = tk.Button(text="Visualize", command=self.VisualizerDo)
        self.VisualizerDo_btn.configure(bg=self.ac1, fg=self.fg)
        self.VisualizerDo_btn.grid(row=15, column=1, pady=(GUI_std_pad, GUI_std_pad))
    
    def chgUniProt_check(self):
        if self.useUniProt.get() == True:
            self.useFeatures1_check.config(state="normal")
            self.useFeatures2_check.config(state="normal")
            self.useFeatures3_check.config(state="normal")
        else:
            self.useFeatures1_check.config(state="disabled")
            self.useFeatures2_check.config(state="disabled")
            self.useFeatures3_check.config(state="disabled")
    
    def set_color_preset(self, arrg):
        global PEPTOGRAPH_COLORS
        global PEPTOGRAPH_COLORS_PRESETS
        del arrg
        PEPTOGRAPH_COLORS = PEPTOGRAPH_COLORS_PRESETS[self.preset_name.get()]*100
        self.DataGroup1_color.configure(bg='#'+PEPTOGRAPH_COLORS[0])
        self.DataGroup2_color.configure(bg='#'+PEPTOGRAPH_COLORS[1])
        self.DataGroup3_color.configure(bg='#'+PEPTOGRAPH_COLORS[2])
        self.DataGroup4_color.configure(bg='#'+PEPTOGRAPH_COLORS[3])
        self.DataGroup5_color.configure(bg='#'+PEPTOGRAPH_COLORS[4])
        self.DataGroup6_color.configure(bg='#'+PEPTOGRAPH_COLORS[5])
        self.DataGroup7_color.configure(bg='#'+PEPTOGRAPH_COLORS[6])
        self.DataGroup8_color.configure(bg='#'+PEPTOGRAPH_COLORS[7])
    
    def chg1(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[0]:
            PEPTOGRAPH_COLORS[0] = this_color[1].lstrip("#")
            self.DataGroup1_color.configure(bg=this_color[1])
    def chg2(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[1] = this_color[1].lstrip("#")
            self.DataGroup2_color.configure(bg=this_color[1])
    def chg3(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[2] = this_color[1].lstrip("#")
            self.DataGroup3_color.configure(bg=this_color[1])
    def chg4(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[3] = this_color[1].lstrip("#")
            self.DataGroup4_color.configure(bg=this_color[1])
    def chg5(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[4] = this_color[1].lstrip("#")
            self.DataGroup5_color.configure(bg=this_color[1])
    def chg6(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[5] = this_color[1].lstrip("#")
            self.DataGroup6_color.configure(bg=this_color[1])
    def chg7(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[6] = this_color[1].lstrip("#")
            self.DataGroup7_color.configure(bg=this_color[1])
    def chg8(self):
        this_color = askcolor(title="Choose color ...")
        if this_color[1]:
            PEPTOGRAPH_COLORS[7] = this_color[1].lstrip("#")
            self.DataGroup8_color.configure(bg=this_color[1])
    
    def add_1(self):
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
        
    def del_1(self):
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
    
    def add_2(self):
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
        
    def del_2(self):
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
    
    def add_3(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup3_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup3_vars.set(value=result_list)
        
    def del_3(self):
        if self.DataGroup3_list.curselection() == ():
            return
        selection = self.DataGroup3_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup3_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup3_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup3_vars.set(value=existing_list)
    
    def add_4(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup4_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup4_vars.set(value=result_list)
        
    def del_4(self):
        if self.DataGroup4_list.curselection() == ():
            return
        selection = self.DataGroup4_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup4_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup4_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup4_vars.set(value=existing_list)
        
    def add_5(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup5_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup5_vars.set(value=result_list)
        
    def del_5(self):
        if self.DataGroup5_list.curselection() == ():
            return
        selection = self.DataGroup5_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup5_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup5_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup5_vars.set(value=existing_list)
    
    def add_6(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup6_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup6_vars.set(value=result_list)
        
    def del_6(self):
        if self.DataGroup6_list.curselection() == ():
            return
        selection = self.DataGroup6_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup6_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup6_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup6_vars.set(value=existing_list)
        
    def add_7(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup7_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup7_vars.set(value=result_list)
        
    def del_7(self):
        if self.DataGroup7_list.curselection() == ():
            return
        selection = self.DataGroup7_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup7_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup7_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup7_vars.set(value=existing_list)
        
    def add_8(self):
        if self.DataSets_list.curselection() == ():
            return
        
        selection = self.DataSets_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataSets_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup8_vars.get()[1:-1].split(',')]

        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list:
                selected_list.remove(i)

        result_list = existing_list + selected_list
        while ("" in result_list): # brez tega mam skos eno prazno polje, bruh
            result_list.remove("")
        
        self.DataGroup8_vars.set(value=result_list)
        
    def del_8(self):
        if self.DataGroup8_list.curselection() == ():
            return
        selection = self.DataGroup8_list.curselection()
        value_list = [line.strip(' \'') for line in self.DataGroup8_vars.get()[1:-1].split(',')]
        selected_list = [value_list[index] for index in selection]
        existing_list = [line.strip(' \'') for line in self.DataGroup8_vars.get()[1:-1].split(',')]
        
        # merge w/o duplicates
        for i in selected_list[:]:
            if i in existing_list: #praktično nepotreben pogoj, ma vseeno
                existing_list.remove(i)
                
        while ("" in existing_list):
            existing_list.remove("")
        self.DataGroup8_vars.set(value=existing_list)
    
    def another_exp_pls(self):
        global resize_factor
        # DataGroup2
        if not self.DataGroup2_list.winfo_ismapped():
            #prvič dobim resize factor
            resize_factor = round(self.winfo_width() / 3) + round(GUI_list_pad / 3)
            
            self.DataGroup2_text.grid(row=1, column=2,sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup2_color.grid(row=1, column=2, sticky='E')
            self.DataGroup2_list.grid(row=2, column=2, padx=(GUI_list_pad,0))
            self.add2.grid(row=3, column=2, padx=(0,GUI_button_pad))
            self.remove2.grid(row=3, column=2, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=3)
            
            self.geometry(str(int(self.winfo_width())+int(self.AddColumn.winfo_width())+2*GUI_list_pad)+"x"+str(self.winfo_height())) #prvič se poveča samo za širino gumba
            return
        # DataGroup3
        if not self.DataGroup3_list.winfo_ismapped():
            self.DataGroup3_text.grid(row=1, column=3, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup3_color.grid(row=1, column=3, sticky='E')
            self.DataGroup3_list.grid(row=2, column=3, padx=(GUI_list_pad,0))
            self.add3.grid(row=3, column=3, padx=(0,GUI_button_pad))
            self.remove3.grid(row=3, column=3, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=4, pady=(10,0))
            
            self.geometry(str(int(self.winfo_width())+resize_factor)+"x"+str(self.winfo_height()))
            return
        # DataGroup4
        if not self.DataGroup4_list.winfo_ismapped():
            self.DataGroup4_text.grid(row=1, column=4, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup4_color.grid(row=1, column=4, sticky='E')
            self.DataGroup4_list.grid(row=2, column=4, padx=(GUI_list_pad,0))
            self.add4.grid(row=3, column=4, padx=(0,GUI_button_pad))
            self.remove4.grid(row=3, column=4, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=5, pady=(10,0))
            
            self.geometry(str(int(self.winfo_width())+resize_factor)+"x"+str(self.winfo_height()))
            return
        # DataGroup5
        if not self.DataGroup5_list.winfo_ismapped():
            self.DataGroup5_text.grid(row=1, column=5, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup5_color.grid(row=1, column=5, sticky='E')
            self.DataGroup5_list.grid(row=2, column=5, padx=(GUI_list_pad,0))
            self.add5.grid(row=3, column=5, padx=(0,GUI_button_pad))
            self.remove5.grid(row=3, column=5, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=6, pady=(10,0))
            
            self.geometry(str(int(self.winfo_width())+resize_factor)+"x"+str(self.winfo_height()))
            return
        # DataGroup6
        if not self.DataGroup6_list.winfo_ismapped():
            self.DataGroup6_text.grid(row=1, column=6, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup6_color.grid(row=1, column=6, sticky='E')
            self.DataGroup6_list.grid(row=2, column=6, padx=(GUI_list_pad,0))
            self.add6.grid(row=3, column=6, padx=(0,GUI_button_pad))
            self.remove6.grid(row=3, column=6, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=7, pady=(10,0))
            
            self.geometry(str(int(self.winfo_width())+resize_factor)+"x"+str(self.winfo_height()))
            return
        # DataGroup7
        if not self.DataGroup7_list.winfo_ismapped():
            self.DataGroup7_text.grid(row=1, column=7, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup7_color.grid(row=1, column=7, sticky='E')
            self.DataGroup7_list.grid(row=2, column=7, padx=(GUI_list_pad,0))
            self.add7.grid(row=3, column=7, padx=(0,GUI_button_pad))
            self.remove7.grid(row=3, column=7, padx=(GUI_button_pad,0))
            
            self.AddColumn.grid_remove()
            self.AddColumn.grid(row=2, column=8, pady=(10,0))
            
            self.geometry(str(int(self.winfo_width())+resize_factor)+"x"+str(self.winfo_height()))
            return
        # DataGroup8
        if not self.DataGroup8_list.winfo_ismapped():
            self.DataGroup8_text.grid(row=1, column=8, sticky='W', padx=(GUI_std_pad,0))
            self.DataGroup8_color.grid(row=1, column=8, sticky='E')
            self.DataGroup8_list.grid(row=2, column=8, padx=(GUI_list_pad,0))
            self.add8.grid(row=3, column=8, padx=(0,GUI_button_pad))
            self.remove8.grid(row=3, column=8, padx=(GUI_button_pad,0))
            self.AddColumn.grid_remove()
            self.geometry(str(int(self.winfo_width())+resize_factor-int(self.AddColumn.winfo_width()))+"x"+str(self.winfo_height()))
            return
    
    def group_normalization_alert(self):
        if self.group_normalization.get() == True:
            tkinter.messagebox.showwarning(title="Group normalization", message="Using this setting is not recommended since intensities will be handeled in each experiment separately!")
        return
    
    def VisualizerDo(self):
        global usr_name
        global DataGroup_names
        global DataGroup_values
        global EXP_N
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
        global data_indexes
        global data_titles
        global setMW_Labels
        global peptide_threshold
        
        EXP_N = 0
        DataGroup_names = []
        DataGroup_values = []
        DEBUG_txt.append(f'NOTICE: Acquiring user settings from GUI...')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        # Experiment name
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
        
        if self.DataGroup1_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup1_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup1_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup2_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup2_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup2_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup3_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup3_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup3_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup4_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup4_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup4_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup5_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup5_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup5_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup6_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup6_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup6_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup7_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup7_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup7_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1
        if self.DataGroup8_vars.get()[1:-1]:
            tmp = ""
            tmp_arr = []
            DataGroup_names.append(self.DataGroup8_name.get())
            tmp = [line.strip(' \'') for line in self.DataGroup8_vars.get()[1:-1].split(',')]
            while ("" in tmp):
                tmp.remove("")
            for t in tmp:
                tmp_arr.append(data_titles.index(t))
            DataGroup_values.append(tmp_arr)
            EXP_N+=1

        #DEBUG
        DEBUG_txt.append(f"NOTICE: EXP_N = {EXP_N}\nDataGroup_names = {json.dumps(DataGroup_names)}\nNOTICE: DataGroup_values = {json.dumps(DataGroup_values)}") #v prejšnjem debugu izpisujem DataGroup2 in ne DataGroup2_names
        if NOTICE:
            print(DEBUG_txt[-1])

        if EXP_N > 0:
            N_max = len(DataGroup_values[0])
            for i in range(EXP_N):
                if len(DataGroup_values[i]) != N_max:
                    print(f"{bcolors.FAIL}ERROR! \"{DataGroup_names[i]}\" is not the same size as \"{DataGroup_names[0]}\"!{bcolors.ENDC}")
                    return
        else:
            print(f"{bcolors.FAIL}ERROR! No data groups defined!{bcolors.ENDC}")
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
        DEBUG_txt.append(f'NOTICE: N_max = {N_max}')
        if NOTICE:
            print(DEBUG_txt[-1])
            
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
        
        setMW = self.setMW.get()
        DEBUG_txt.append(f'NOTICE: setMW = {setMW}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        group_normalization = self.group_normalization.get()
        DEBUG_txt.append(f'NOTICE: group_normalization = {group_normalization}')
        if NOTICE:
            print(DEBUG_txt[-1])
        
        self.quit() #je to sploh potrebno, glede na to, da uničim na koncu
        self.destroy()

# Vnos molekulskih mass
class PeptideVisualizerSetMwGUI(tk.Tk):
    def __init__(self, *args, **kwargs):
        global N_max
        global EXP_N
        protein_length_example = 255
        tk.Tk.__init__(self, *args, **kwargs)
        
        # DARK/LIGHT Mode
        if DARK_MODE:
            self.bg, self.fg, self.ac1, self.ac2 = ('#282828', 'white', '#404040', '#B3B3B3')
        else:
            self.bg, self.fg, self.ac1, self.ac2 = ('#f0f0ed', 'black', '#ffffff', '#E9DAC1')
        wh = round(900 * OS_scale)
        ww = round(550 * OS_scale)
       
        self.tk.call('tk', 'scaling', scale_f)
        #self.font = ('Serif', round(12*scale_f))
        self.configure(bg=self.bg, borderwidth=0)
        self.title('Set Y-ax labels')
        self.geometry(f"{wh}x{ww}")
        self.resizable()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [3, 1]})
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
        
        plt.subplots_adjust(left = 0.1, bottom = 0.11, right = 0.96, top = 0.87, wspace=0)
        
        canvas = FigureCanvasTkAgg(fig, master = self)  
        canvas.draw()
        canvas.get_tk_widget().place(x=0, y=0)
        
        self.setMWlabels = list()
        self.setMWtexts = list()
        spacing_factor = 355 / N_max
        
        for i in range(N_max-1):
            self.setMWtext = tk.Entry(width=5)
            self.setMWtext.place(x=round(50*OS_scale), y=round((55+((i+1)*spacing_factor))*OS_scale)) # to včasih zataji. Vsaj na applu
            #self.setMWtext.place(x=50, y=(90+(spacing_factor/2)+(spacing_factor*i)))
            self.setMWtexts.append(self.setMWtext)
        
        #presets
        MW_options = []
        for preset_name in MW_PRESET.keys():
            if N_max in MW_PRESET[preset_name]:
                MW_options.append(preset_name)
        if len(MW_options) != 0:
            self.preset_name = tk.StringVar(value="Gel presets")
            self.preset_menu = tk.OptionMenu(self, self.preset_name, *MW_options, command=self.set_preset)
            self.preset_menu.config(bg=self.bg, fg=self.fg, activebackground=self.ac1, activeforeground=self.fg, borderwidth=0)
            self.preset_menu["menu"].config(bg=self.bg, fg=self.fg)
            self.preset_menu.place(x=round(100*OS_scale), y=round(507*OS_scale))
        
        self.setYax_btn = tk.Button(text="Submit", command=self.setYax)
        self.setYax_btn.configure(bg=self.ac1, fg=self.fg)
        self.setYax_btn.place(x=round(420*OS_scale), y=round(510*OS_scale))
    
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
    # EXP_N
    file.write("EXP_N: ")
    try:
        EXP_N
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(str(EXP_N)+"\n")
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
    # DataGroup_names
    file.write("DataGroup_names: ")
    try:
        DataGroup_names
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(DataGroup_names)+"\n")
    # DataGroup_values
    file.write("DataGroup_values: ")
    try:
        DataGroup_values
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(DataGroup_values)+"\n")
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
    # Total_intensities
    file.write("Total_intensities: ")
    try:
        Total_intensities
    except NameError:
        file.write("NOT DEFINED!\n")
    else:
        file.write(json.dumps(Total_intensities)+"\n")
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
    # konec
    file.close()

##############################################################################
#
#                       0. ONLINE CHECK + UPDATE MODUL
#
##############################################################################

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
    root = tk.Tk()
    root.withdraw()
    root.update()
    peptides_file = askopenfilename(title = "Select peptides file",filetypes=[("peptides.txt", ".txt .csv")])
    root.destroy()
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
        
        # try: ?!
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

# GUI jupiii ToDO: CLI rabi prnovo!
if GUI:
    App = PeptideVisualizerGUI()
    App.mainloop()
    #App.destroy()
    if not DataGroup_values:
        print(f"{bcolors.FAIL}ERROR! Experiment settings not set! Terminating ...{bcolors.ENDC}")
        pause()
        exit()
    if setMW and N_max > 1:
        App_setMW = PeptideVisualizerSetMwGUI()
        App_setMW.mainloop()
        #App_setMW.destroy()
        if setMW_Labels == ['']*(N_max-1):
            setMW = False
            DEBUG_txt.append(f"WARNING! No Y-ax labels given. Y-ax will be labeled by slice number. setMW = {setMW}")
            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
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

peptide_threshold *= 0.01

if SVG_OUTPUT:
    OUTPUT_EXTENTION = '.svg'
else:
    OUTPUT_EXTENTION = '.png'

##############################################################################
#               peptides.txt handling
##############################################################################
times = {}
times[0] = time.time()
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
    
times[1] = time.time()
DEBUG_txt.append(f"{bcolors.OKBLUE}Read {line_count} peptides in {round((times[1] - times[0])*1000)} ms{bcolors.ENDC}")
print(DEBUG_txt[-1])

#prepare list for processing
multi_IDs = []
for id in DB:
    multi_IDs.append(id) 

DEBUG_txt.append(f'All the required data is stored in the memory. Starting peptide visualization of {len(DB.keys())} proteins ...')
if NOTICE:
    print(DEBUG_txt[-1])

##############################################################################
#               FASTA handling
##############################################################################
if useUniProt:
    DEBUG_txt.append(f" ( 1/3 ) Retriving FASTA from UniProt ...")
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
    DEBUG_txt.append(f" ( 1/3 ) Reading FASTA ...")
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

times[2] = time.time()
DEBUG_txt.append(f"{bcolors.OKBLUE}handled {len(FASTA)} FASTA sequences in {round(times[2] - times[1])} s{bcolors.ENDC}")
print(DEBUG_txt[-1])

##############################################################################
#               Statistična obdelava rezultatov
##############################################################################

DEBUG_txt.append(f" ( 2/3 ) Statistical analysis ... ")
print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")

# Deklariram spremenljivke, da se ne deklarirajo v zankah
NOZMALIZATION_FACTORS = {}
PEPTIDES_quantitative_avg = {}
PEPTIDES_quantitative_err = {}
PEPTIDES_quantitative_avg_norm = {}
SLICE_quantitative = {}
SLICE_quantitative_err = {}
Total_intensities = {}
K_factor = {}
SUMMARY = {}
#for this_protein, these_peptides in DB.items():
    #SLICE_quantitative[this_protein] = {}
    #SLICE_quantitative_err[this_protein] = {}
    #Total_intensities[this_protein] = {}
    #NOZMALIZATION_FACTORS[this_protein] = {}
    #for exp_ID in range(EXP_N):
        #SLICE_quantitative[this_protein][exp_ID] = {}
        #SLICE_quantitative_err[this_protein][exp_ID] = {}
        #Total_intensities[this_protein][exp_ID] = 0
    #for this_peptide in these_peptides:
        #PEPTIDES_quantitative_avg[this_peptide] = {}
        #PEPTIDES_quantitative_err[this_peptide] = {}
        #PEPTIDES_quantitative_avg_norm[this_peptide] = {}
        #for exp_ID in range(EXP_N):
            #PEPTIDES_quantitative_avg[this_peptide][exp_ID] = {}
            #PEPTIDES_quantitative_err[this_peptide][exp_ID] = {}
            #PEPTIDES_quantitative_avg_norm[this_peptide][exp_ID] = {}

# Dejanska statistična obdelava
with alive_bar(len(DB.keys()), dual_line=True, theme='classic') as bar:
    for this_protein, these_peptides in DB.items():
        bar.text = f'Statistically analyzing protein {this_protein}'
        
        SLICE_quantitative[this_protein] = {}
        SLICE_quantitative_err[this_protein] = {}
        Total_intensities[this_protein] = {}
        NOZMALIZATION_FACTORS[this_protein] = {}
        
        SUMMARY[this_protein] = [slugify(this_protein.replace("|", "_")), 0, 0]
        # SUMMARY[protein_id][0] = Protein name
        # SUMMARY[protein_id][1] = sum(intensity)
        # SUMMARY[protein_id][2] = sum(peptides)
        # SUMMARY[protein_id][3] = coverage - dodan kasneje
        # SUMMARY[protein_id][4] = K-factor - dodan kasneje
        
        # FASTA check
        try: # to ne dela, ker se sesuje. morda bi lahko prečekiral pri fasti, če je kompletna?
            FASTA[this_protein]
        except NameError:
            print(f'{bcolors.FAIL}ERROR! FASTA entry for "{this_protein}" does not exist! Check FASTA file. Terminating ...{bcolors.ENDC}')
            pause()
            exit(2)

        for exp_ID in range(EXP_N): # zapeljem čez vse eksperimente
            SLICE_quantitative[this_protein][exp_ID] = {}
            SLICE_quantitative_err[this_protein][exp_ID] = {}
            Total_intensities[this_protein][exp_ID] = 0
            NOZMALIZATION_FACTORS[this_protein][exp_ID] = []
            ###
            # SLICE_quantitative in SLICE_quantitative_err, dokler še imam intenzitete absolutne.
            ###
            slice = 1
            for i in range(N_max): # zapeljem čez vse rezine
                #ID = DataGroup_values[exp_ID][i]
                replicate_intensity = 0
                work_list = []
                for j in range(replicates):
                    rep_ID = DataGroup_values[exp_ID][i+(j*N_max)]
                    for this_peptide in these_peptides:
                        replicate_intensity += PEPTIDES_quantitative[this_peptide][rep_ID]
                    work_list.append(replicate_intensity)
    
                if replicates > 1:
                    SLICE_quantitative[this_protein][exp_ID][slice] = statistics.mean(work_list)
                    SLICE_quantitative_err[this_protein][exp_ID][slice] = statistics.stdev(work_list)
                else:
                    SLICE_quantitative[this_protein][exp_ID][slice] = work_list[0]
                    SLICE_quantitative_err[this_protein][exp_ID][slice] = 0
                del work_list
                del replicate_intensity
                slice += 1
                #Totalne intenzitete so pripisane pred normalizacijo rezin
                Total_intensities[this_protein][exp_ID] += SLICE_quantitative[this_protein][exp_ID][i+1]
        
            ###
            # transformiram PEPTIDES_quantitative v PEPTIDES_quantitative_avg in PEPTIDES_quantitative_err pa pol uporabljam ta dva kasneje bojo normalizirani
            ###
            for this_peptide in these_peptides:
                try: # daj odstrani to iz loopa
                    PEPTIDES_quantitative_avg[this_peptide]
                except KeyError:
                    PEPTIDES_quantitative_avg[this_peptide] = {}
                    PEPTIDES_quantitative_err[this_peptide] = {}
                    PEPTIDES_quantitative_avg_norm[this_peptide] = {}
                
                try:
                    PEPTIDES_quantitative_avg[this_peptide][exp_ID]
                except KeyError:
                    PEPTIDES_quantitative_avg[this_peptide][exp_ID] = {}
                    PEPTIDES_quantitative_err[this_peptide][exp_ID] = {}
                    PEPTIDES_quantitative_avg_norm[this_peptide][exp_ID] = {}
                
                if replicates > 1:
                    avg_list = []
                    err_list = []
                    for i in range(N_max): #grem čez rezine
                        ID = DataGroup_values[exp_ID][i]
                        work_list = []
                        for j in range(replicates):
                            rep_ID = DataGroup_values[exp_ID][i+(j*N_max)]
                            work_list.append(PEPTIDES_quantitative[this_peptide][rep_ID])
                        
                        avg_list.append(statistics.mean(work_list))
                        if avg_list[i] > 0: # da ne pride do deljenja z 0
                            err_list.append(statistics.stdev(work_list))
                        else:
                            err_list.append(0)
                        del(work_list)

                    PEPTIDES_quantitative_avg[this_peptide][exp_ID] = avg_list
                    PEPTIDES_quantitative_err[this_peptide][exp_ID] = err_list
                    del avg_list
                    del err_list
                
                else: # če ni replikatov je avg samo prepisana vrednost
                    for i in range(N_max):
                        ID = DataGroup_values[exp_ID][i]
                        PEPTIDES_quantitative_avg[this_peptide][exp_ID][i] = PEPTIDES_quantitative[this_peptide][ID]
                    PEPTIDES_quantitative_err[this_peptide][exp_ID] = [0] * N_max
                
            ###
            # Normalizacija peptidov (PEPTIDES_quantitative_avg) + SUMMARY
            ###
            if group_normalization: # Fukjena normalizacija
                tmp_list=[]
                for this_peptide in these_peptides:
                    SUMMARY[this_protein][2] += 1 #preštejem peptide
                    for i in range(N_max):
                        tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][exp_ID][i])
                #SUMMARY
                SUMMARY[this_protein][1] += sum(tmp_list) #seštejem intenziteto

                normalization_factor = max(tmp_list)
                if normalization_factor > 0:
                    NOZMALIZATION_FACTORS[this_protein][exp_ID] = normalization_factor
                    for this_peptide in these_peptides:
                        work_list = []
                        for i in range(N_max):
                            work_list.append(PEPTIDES_quantitative_avg[this_peptide][exp_ID][i] / normalization_factor)
                        PEPTIDES_quantitative_avg_norm[this_peptide][exp_ID][i] = work_list
                        del work_list
                else:
                    NOZMALIZATION_FACTORS[this_protein][exp_ID] = False
                    for this_peptide in these_peptides:
                        PEPTIDES_quantitative_avg_norm[this_peptide] = [0] * N_max
                del normalization_factor
                del tmp_list
        
        if not group_normalization: # Normalna normalizacija
            tmp_list=[]
            for this_peptide in these_peptides:
                SUMMARY[this_protein][2] += 1 #preštejem peptide
                for exp_ID in range(EXP_N):
                    for i in range(N_max):
                        tmp_list.append(PEPTIDES_quantitative_avg[this_peptide][exp_ID][i])
            #SUMMARY
            SUMMARY[this_protein][1] += sum(tmp_list) #seštejem intenziteto
            
            normalization_factor = max(tmp_list)
            if normalization_factor > 0:
                NOZMALIZATION_FACTORS[this_protein][exp_ID] = normalization_factor #kljub temu, da je za vse eksperimente enak, shranim vsakega posebaj, zaradi tega da bo polje konsistentno s fukjeno normalizacijo
                for this_peptide in these_peptides:
                    for exp_ID in range(EXP_N):
                        for i in range(N_max):
                            PEPTIDES_quantitative_avg_norm[this_peptide][exp_ID][i] = PEPTIDES_quantitative_avg[this_peptide][exp_ID][i] / normalization_factor
            else: #če je normalčizacijski faktor enak 0 pomeni da so vsi elementi enaki 0
                NOZMALIZATION_FACTORS[this_protein] = [False] * EXP_N
                for this_peptide in these_peptides:
                    for exp_ID in range(EXP_N):
                        for i in range(N_max):
                            PEPTIDES_quantitative_avg_norm[this_peptide][exp_ID][i] = 0
            del normalization_factor
            del tmp_list
            
            
        # normalizacija rezin (intenzitete rezine)
        tmp_list = []
        for exp_ID in range(EXP_N):
            tmp_list += SLICE_quantitative[this_protein][exp_ID].values()
        slice_normalization_factor = max(tmp_list)
        if slice_normalization_factor > 0:
            for exp_ID in range(EXP_N):
                for n in range(len(SLICE_quantitative[this_protein][exp_ID].keys())):
                    SLICE_quantitative[this_protein][exp_ID][n+1] /= slice_normalization_factor
                    SLICE_quantitative_err[this_protein][exp_ID][n+1] /= slice_normalization_factor
        del slice_normalization_factor
        
        # Normalizacija totalnih eksperimentov (totalne intenzitete)
        total_normalization_factor = max(Total_intensities[this_protein].values())    
        if total_normalization_factor != 0:
            for exp_ID in range(EXP_N):
                Total_intensities[this_protein][exp_ID] /= total_normalization_factor
        del total_normalization_factor
        
        bar() #progress
        # DEBUG
        #print(f"\n\nPEPTIDES_quantitative_avg = {PEPTIDES_quantitative_avg} \n\nPEPTIDES_quantitative_avg_norm = {PEPTIDES_quantitative_avg_norm} \n\nNOZMALIZATION_FACTORS = {NOZMALIZATION_FACTORS} \n\nPEPTIDES_quantitative_err = {PEPTIDES_quantitative_err} \n\nSLICE_quantitative = {SLICE_quantitative} \n\nSLICE_quantitative_err = {SLICE_quantitative_err} \n\nTotal_intensities = {Total_intensities}")
        #pause()

#write_debug()
#pause()

times[3] = time.time()
DEBUG_txt.append(f"{bcolors.OKBLUE}Statistically analysed {len(FASTA)} proteins in {round((times[3] - times[2])*1000)} ms{bcolors.ENDC}")
print(DEBUG_txt[-1])

##############################################################################
#               Generiranje rezultatov
##############################################################################

DEBUG_txt.append(f' ( 3/3 ) Generating PROTOMAP results ...')
print(f"\n{bcolors.BOLD}{DEBUG_txt[-1]}{bcolors.ENDC}")

if not os.path.isdir(work_dir):
    os.mkdir(work_dir) # ToDo: test za pisanje? To bi v bistvu moralo bit že prej... Also premakni generiranje HTML pred to srce

###
# Risanje cbar.svg
###
def drawColorbars():
    global EXP_N
    
    PEPTOGRAPH_COLORS_MATRIX = {}
    PEPTOGRAPH_COLORMAP = {}
    PEPTOGRAPH_sm = {}
    PEPTOGRAPH_Cbar = {}
    # Barve za cmap
    N = 256
    for i in range(EXP_N):
        PEPTOGRAPH_COLORS_MATRIX[i] = np.ones((N, 4))
        for j in range(3):
            j2 = j*2
            dec_color = int(PEPTOGRAPH_COLORS[i][j2:j2+2], 16)
            PEPTOGRAPH_COLORS_MATRIX[i][:, j] = np.linspace(dec_color/N, dec_color/N, N)
            del dec_color
        PEPTOGRAPH_COLORS_MATRIX[i][:, 3] = np.linspace(alpha_boost, 1, N)
        PEPTOGRAPH_COLORMAP[i] = ListedColormap(PEPTOGRAPH_COLORS_MATRIX[i])
    # Padding fix za majhne EXP_N
    size_x = 0.5+EXP_N/3.7
    pad_r = 0.5
    if EXP_N > 1 and EXP_N < 8:
        pad_r+=0.3*EXP_N/8
    else:
        pad_r = 0.8
        
    #Narišem
    cBarPlt, ax = plt.subplots(1,EXP_N,figsize=(size_x,5))
    norm = plt.Normalize(alpha_boost, 1)
    
    for i in range(EXP_N):
        PEPTOGRAPH_sm[i] = plt.cm.ScalarMappable(cmap=PEPTOGRAPH_COLORMAP[(EXP_N-1-i)], norm=norm) #ToDo: kaj če je i večji kot je bav? naj to nekje preverja oz. poveča matrix (EXP_N-1-i) zato, da bo čisto na levi prvi exp
        PEPTOGRAPH_sm[i].set_array([])
        
        if EXP_N > 1:
            PEPTOGRAPH_Cbar[i] = plt.colorbar(PEPTOGRAPH_sm[i], ax=ax[(EXP_N-1-i)], ticks=[alpha_boost, 1], aspect=22, fraction=1)
            ax[i].set_axis_off()
        else:
            PEPTOGRAPH_Cbar[i] = plt.colorbar(PEPTOGRAPH_sm[i], ax=ax, ticks=[alpha_boost, 1], aspect=22, fraction=1)
            ax.set_axis_off()
        PEPTOGRAPH_Cbar[i].ax.text(0.52, 0.1, DataGroup_names[(EXP_N-1-i)], rotation='vertical', va='baseline', ha='center')
        
        if i == 0:
            PEPTOGRAPH_Cbar[i].set_ticklabels([str(round(alpha_boost*100))+"%", str(100)+"%"])
        else:
            PEPTOGRAPH_Cbar[i].set_ticks([])
        
    cBarPlt.subplots_adjust(left = 0.05, bottom = 0.11, right = pad_r, top = 0.87, wspace=0)
    cBarPlt.savefig(work_dir+'cbar'+OUTPUT_EXTENTION)
    return True

###
# Risanje legend.svg
###
def drawLegend():
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
    return True


# Dejansko narišem vse
drawColorbars()
if ONLINE and useUniProt_features:
    drawLegend()

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

peptograph_exp_spacing /= 100
peptograph_exp_height = (1 - ((EXP_N + 1)*peptograph_exp_spacing)) / EXP_N
def drawPeptographs(protein):
    #SUMMARY
    total_coverage = [False] * len(FASTA[protein]) #ToDo: že prej preverjanje če FASTA[protein] sploh obstaja
    #PRIPRAVIM FAJLE
    result_file = work_dir+slugify(protein.replace("|", "_"))+'.txt' #TXT
    result_html = work_dir+slugify(protein.replace("|", "_"))+'.html' #HTML
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
            print(f"{bcolors.FAIL}ERROR! Given fasta file doesn't match fasta file provided to MaxQuant!\nERROR found at entry \"{protein}\"!{bcolors.ENDC}") #ToDo: naj ponudi za to rešitev
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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [3, 1]})
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
        
        # narišem skupno intenziteto v eksperimentu
        if EXP_N > 1: # če je samo ena skupina je res pointles in čudno zgleda en bervni kvadrat
            total_height = (uniProtFeature_coeficient*(N_max+1) - peptograph_exp_spacing*(EXP_N+1)) / EXP_N
            for exp_ID in range(EXP_N):
                this_color = "#"+PEPTOGRAPH_COLORS[exp_ID]
                y_pos = N_max+0.5+total_height/2+peptograph_exp_spacing+exp_ID*(total_height+peptograph_exp_spacing)
                ax2.barh(y_pos, Total_intensities[protein][exp_ID], height=total_height, color=this_color, alpha=0.5)
                #print(f"Rišem skupno intenziteto pri y_pos={y_pos}")
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
        l_offest = 0.1 + (label_length_delta * 0.009)
    else:
        l_offest = 0.1
        
    l_offest -= round(EXP_N/100)
    r_offset = 0.96 + round(EXP_N/200)
    plt.subplots_adjust(left = l_offest, bottom = 0.11, right = r_offset, top = 0.87, wspace=0)
    
    # Združen coverage
    COVERAGE_EXP = {} #da bi bolj pregledno
    for exp_ID in range(EXP_N):
        COVERAGE_EXP[exp_ID] = list(empty_char * len(FASTA[protein]))
    
    #Začnem obdelavo posameznih slajsov in eksperimentov
    slice = 1
    for i in range(N_max):
        #DEBUG_txt.append(f'NOTICE: working on slice {slice} in protein {protein} ...')
        COVERAGE_LIST = {}
        COVERAGE = {}

        for exp_ID in range(EXP_N):
            COVERAGE_LIST[exp_ID] = list(empty_char * len(FASTA[protein]))
            this_color = "#"+PEPTOGRAPH_COLORS[exp_ID]
            
            for peptide in DB[protein]:
                peptide_start = PEPTIDES_qualitative[peptide][0]
                peptide_length = PEPTIDES_qualitative[peptide][1]
                
                try:
                    peptide_intensity = PEPTIDES_quantitative_avg_norm[peptide][exp_ID][i]
                except:
                    DEBUG_txt.append(f'ERROR: Peptide "{peptide}" not found in variable PEPTIDES_quantitative_avg_norm[{peptide}][{exp_ID}][{i}]!')
                    write_debug()
                    print(bcolors.FAIL+DEBUG_txt[-1]+bcolors.ENDC)
                    pause()
                    exit()
                 
                alpha_intensity = (1-alpha_boost)*peptide_intensity + alpha_boost
                y_pos = slice-0.5+peptograph_exp_spacing+(peptograph_exp_spacing*exp_ID)+(peptograph_exp_height*exp_ID)
                #print(f"Rišem {this_peptide} i={peptide_intensity} y_pos={y_pos}")
                if peptide_intensity > peptide_threshold:
                    # GRAPH(SVG) OUTPUT MAIN
                    ax1.add_patch(Rectangle((peptide_start, y_pos), peptide_length, peptograph_exp_height, alpha=alpha_intensity, linewidth=0, color=this_color))
                    # TXT RESULTS
                    if TXT_OUTPUT:
                        if peptide_start + peptide_length <= len(COVERAGE_LIST[exp_ID]):
                            for n in range(peptide_length):
                                total_coverage[peptide_start + n] = True
                                if n > 0 and n < (peptide_length-1):
                                    COVERAGE_LIST[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_LIST[exp_ID][peptide_start + n], 'mid')
                                    COVERAGE_EXP[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_EXP[exp_ID][peptide_start + n], 'mid')
                                elif n == 0:
                                    COVERAGE_LIST[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_LIST[exp_ID][peptide_start + n], 'start')
                                    COVERAGE_EXP[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_EXP[exp_ID][peptide_start + n], 'start')
                                elif n == (peptide_length-1):
                                    COVERAGE_LIST[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_LIST[exp_ID][peptide_start + n], 'end')
                                    COVERAGE_EXP[exp_ID][peptide_start + n] = sequence_adjust(COVERAGE_EXP[exp_ID][peptide_start + n], 'end')
                        else:
                            DEBUG_txt.append(f"WARNING: peptide \"{peptide}\" does not fit onto protein \"{protein}\"! Disregarding it in \"{data_titles[DataGroup_values[exp_ID][i]]}\" ({PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1]} > {len(COVERAGE_LIST[exp_ID])})")
                            print(f"{bcolors.WARNING}{DEBUG_txt[-1]}{bcolors.ENDC}")
                elif peptide_intensity > 0 and NOTICE:
                    #DEBUG_txt.append(f"NOTICE: peptide \"{peptide}\" falls out of {round(peptide_threshold*100, 0)}% threshold ({round(peptide_intensity*100, 3)}%). Disregarding it in {protein} {data_titles[DataGroup_values[exp_ID][i]]}")
                    #print(DEBUG_txt[-1])
                    pass #enostavno je preveč outputa. Morda če bi mel nek verbos mode
            
            # dodam posamezni slice in experiment v txt in narišem skupno intenziteto slica + napako, če je podana
            if TXT_OUTPUT:
                COVERAGE[exp_ID] = ''.join(COVERAGE_LIST[exp_ID])
                result.write(data_titles[DataGroup_values[exp_ID][i]]+"\t"+COVERAGE[exp_ID]+"\n")
                del COVERAGE_LIST[exp_ID]
        
            # rišem dve intenziteti
            if error_bar and replicates > 1:
                this_err = SLICE_quantitative_err[protein][exp_ID][slice]
                err_kwargs = dict(elinewidth=err_width/EXP_N, capsize=err_cap/EXP_N, capthick=err_capthick/EXP_N)
            else:
                this_err = 0
                err_kwargs = dict(elinewidth=0, capsize=0, capthick=0)
              
            ax2.barh(y_pos+peptograph_exp_height/2, SLICE_quantitative[protein][exp_ID][slice], xerr=this_err, height=peptograph_exp_height, color=this_color, alpha=0.5, ecolor='black', error_kw=err_kwargs)

        slice+=1
      
    #SUMMARY
    SUMMARY[protein].append(round((sum(total_coverage) / len(total_coverage))*100, 1)) #total coverage
    
    #ToDo!
    #if DataGroup2:
    #    SUMMARY[protein].append(round(K_factor[protein]*(sum(total_coverage) / len(total_coverage)), 2)) # K-factor - morda ga bi lahko množil s pokritostjo, da ne bo polno proteinov, ki so samo malo zaznani?
    #else:
    #    SUMMARY[protein].append(SUMMARY[protein][-1]) # če ni dataset2, potem izpisuje coverage
    SUMMARY[protein].append(SUMMARY[protein][-1])
    
    #TXT KONEC
    if TXT_OUTPUT:
        result.close()
        # Nardim še .html iz .txt
        #file = open(result_file, 'r')
        #seq = file.readline().split('\t')[1]
        #coverages = []
        #for i in range(N_max):
        #    for exp_ID in range(EXP_N): # še en loop, ker je pokritost zdaj dejansko narjena
        #        coverages.append(file.readline().rstrip('\n').split('\t')[1])
        #file.close()
        
        seq = FASTA[protein]
        for exp_ID in range(EXP_N):
            COVERAGE_EXP[exp_ID] = ''.join(COVERAGE_EXP[exp_ID])
        result = codecs.open(result_html, "w", "utf-8")
        result.write("<table style='font-family: monospace; white-space:pre; margin: auto;'>")
        segment = 0
        AA_n = 0
        while segment*AA_PER_LINE < len(FASTA[protein]):
            numbering_line = ""
            while AA_n < (segment+1)*AA_PER_LINE:
                numbering_line += str(AA_n+1)+" "*(AA_STEP-len(str(AA_n+1)))
                AA_n += AA_STEP
            result.write(f"<tr><td>#</td><td>{numbering_line}</td></tr>")
            result.write(f"<tr><td>SEQ</td><td>{seq[segment*AA_PER_LINE:(segment+1)*AA_PER_LINE]}<td></tr>")
            
            # Eksperimenti, slici so združeni
            for exp_ID in range(EXP_N):
                result.write(f"<tr><td>{DataGroup_names[exp_ID]}</td><td>{COVERAGE_EXP[exp_ID][segment*AA_PER_LINE:(segment+1)*AA_PER_LINE]}<td></tr>")
            # Vsak eksperiment in vsak slice posebej
            #j=0
            #for i in range(N_max):
            #    for exp_ID in range(EXP_N):
            #        result.write(f"<tr><td>{data_titles[DataGroup_values[exp_ID][i]]}</td><td>{coverages[j][segment*AA_PER_LINE:(segment+1)*AA_PER_LINE]}<td></tr>")
            #        j+=1
            segment+=1
            
        result.write("</table>")
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

times[4] = time.time()
DEBUG_txt.append(f"{bcolors.OKBLUE}Drawn {len(FASTA)} peptographs in {round((times[3] - times[2])/60,1)} min{bcolors.ENDC}")
print(DEBUG_txt[-1])

###
##### Generiranje HTML datoteke z rezultati
###
DEBUG_txt.append(f'ALL DONE! Generating {exp_name}.html results file in "{os.path.dirname(peptides_file)}" ... ')
print("\n"+DEBUG_txt[-1])
filenames = []
#for protein,summary_data in SUMMARY.items():
for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][1], reverse=True):
    filenames.append(summary_data[0])

groups_summary = ""
for exp_ID in range(EXP_N):
    groups_summary += DataGroup_names[exp_ID] + ": "
    for data_key in DataGroup_values[exp_ID]:
        groups_summary += data_titles[data_key]+", "
    groups_summary = groups_summary[0:-2]+"<br>"
groups_summary += "Group normalization? "+str(group_normalization)

HTML_CODE = '''
<!DOCTYPE html>
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
		<script src="https://code.jquery.com/jquery-3.6.1.min.js"></script>
		<script>
			const PROTEINS = '''+json.dumps(filenames)+''';
			N = PROTEINS.length - 1;
			last_call = -1;
			
			function loadProtein(ID) {
				last_call = ID;
				txt_file = "'''+exp_name+'''/"+PROTEINS[ID]+".txt";
				html_file = "'''+exp_name+'''/"+PROTEINS[ID]+".html";
				img_file = "'''+exp_name+'''/"+PROTEINS[ID]+"'''+OUTPUT_EXTENTION+'''";
				cbar_file = "'''+exp_name+"/cbar"+OUTPUT_EXTENTION+'''";
				legend_file = "'''+exp_name+'''/legend'''+OUTPUT_EXTENTION+'''";
				//alert(txt_file+"\\n"+html_file+"\\n"+img_file+"\\n"+cbar_file+"\\n"+legend_file)
				document.getElementById("content").innerHTML = '<img src='+img_file+' alt="'+PROTEINS[ID]+'" onClick=""><img src='+cbar_file+' alt="Colorbars"><br><br><object width="100%" height="80%" data="'+html_file+'"></object><br><br><a href="https://www.uniprot.org/help/sequence_annotation"><img src='+legend_file+'></a>';
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
						<th class="main_colum">Protein Name</th>
						<th class="sub_colum">K-factor</th>
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
        HTML_CODE += '''					<tr onclick="loadProtein('''+str(n)+''')">
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
        HTML_CODE += '''					<tr onclick="loadProtein('''+str(n)+''')">
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
os.startfile(os.path.dirname(peptides_file)+'/'+exp_name+".html")

###
##### Konec + DEBUG
###
times[5] = time.time()
DEBUG_txt.append(f"DONE!")
print(f"\n{bcolors.OKGREEN}{DEBUG_txt[-1]}{bcolors.ENDC}")
DEBUG_txt.append(f"Script successfully finished in {round((times[5]-times[0])/60,1)} min.")
print(f"{bcolors.OKBLUE}{DEBUG_txt[-1]}{bcolors.ENDC}\n")

if DEBUG:
    write_debug()
pause()