# (c) 2022
#     Matej Kolarič (matej.kolaric@ijs.si)
#     Biochemistry, molecular and structural biology - B1, Jozef Stefan Institute, Jamova cesta 39, Ljubljana Slovenia
#     under GPLv3 (http://www.gnu.org/licenses/gpl-3.0.html)

##############################################################################
#               Nastavitve
##############################################################################

DEBUG = False
NOTICE = False
bar_height = 0.9
bar_width = 0.8
alpha_boost = 0.05
default_peptide_threshold = 2 # v odstotkih
replicates_MAX = 10000

empty_char = ' '
start_char = '├'
mid_char = '-'
end_char = '┤'
cross_char = '┼'
PROGRAM_NAME = 'PeptideVisulizer'

##############################################################################
#               Moduli
##############################################################################

print("Initilizing "+PROGRAM_NAME+" ...")

import os
os.system("") #colorfix
import csv
import subprocess
import sys
import re
import time
from datetime import datetime
from tkinter.filedialog import askopenfilename
from tkinter import messagebox
import unicodedata
import codecs
import json

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

def pause():
    programPause = input("Press the <ENTER> key to continue...")

# REQUIRED: numpy
try:
    import numpy as np
except ImportError:
    print(bcolors.WARNING+'WARNING: Python module numpy not found! PeptideVisulizer will try to install it for you. If the script fails run "pip install numpy" in command promt.'+bcolors.ENDC)
    pause()
    print("Istalling numpy ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'numpy', "--quiet"])
    time.sleep(1)
    import numpy as np
    print("PeptideVisulizer successfully installed matplotlib v"+numpy.__version__+"\n")

# REQUIRED: matplotlib
try:
    import matplotlib
except ImportError:
    print(bcolors.WARNING+'WARNING: Python module matplotlib not found! PeptideVisulizer will try to install it for you. If the script fails run "pip install matplotlib" in command promt.'+bcolors.ENDC)
    pause()
    print("Istalling matplotlib ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'matplotlib', "--quiet"])
    time.sleep(1)
    import matplotlib
    print("PeptideVisulizer successfully installed matplotlib v"+matplotlib.__version__+"\n")
    
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
matplotlib.use("Agg") #brez te vrstice mi odleti vn iz memory limita (2GB). Agg, is a non-interactive backend that can only write to files.
# RECOMMENDED: alive_progress       ToDo: make it optional
try:
    from alive_progress import alive_bar
except ImportError:
    print(bcolors.WARNING+"WARNING: Python module alive-progress not found! PeptideVisulizer will try to install it for you. If the script fails run \"pip install alive-progress\" in command promt."+bcolors.ENDC)
    pause()
    print("Istalling alive-progress ... ")
    subprocess.check_call([sys.executable, "-m", "pip", "install", 'alive-progress', "--quiet"])
    time.sleep(1)
    from alive_progress import alive_bar
    print("PeptideVisulizer successfully installed alive_progress\n")

##############################################################################
#               Funkcije
##############################################################################

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

def userGrouping(data):    
    DataGroup1 = []
    DataGroup2 = []
    print(bcolors.OKBLUE+"\t*** Detected experiments: ***\nID\t=>NAME\n"+bcolors.ENDC)
    for exp_key in range(len(data)):
        print (exp_key+1, "\t=> ", data[exp_key])

    pattern = re.compile('(\d+)\-(\d+)') # za spreminjanje 1-5 v 1,2,3,4,5
    ##### DataGroup1 #####
    usr_in = input("\nDataGroup1 CSV: ").replace(" ", "").replace("\t", "").replace("\n", "")
    if not usr_in:
        print(bcolors.FAIL+'ERROR! No user input'+bcolors.ENDC)
        pause()
        exit(2)

    val = re.sub(pattern, replace_range_by_ints, usr_in).split(",") # ločim vejice in spremenim - v seznam
    for single in val:
        try:
            data_id = int(single)-1
        except ValueError:
            print(bcolors.FAIL+"ERROR! Not a valid input!"+bcolors.ENDC)
            pause()
            exit(2)
        
        if data_id >= 0 and data_id <= len(data):
            if data_id not in DataGroup1:
                DataGroup1.append(data_id)
            else:
                print(bcolors.WARNING+"WARNING: '"+str(data_id+1)+"' is entered twice! Removing second instance!"+bcolors.ENDC)
        else:
            print(bcolors.WARNING+"WARNING: '"+str(data_id+1)+"' is out of range! Removing it!"+bcolors.ENDC)

    #print(DataGroup1)
    print(bcolors.OKBLUE+"DataGroup1: ", end = '')
    for key in DataGroup1:
        print(data[key]+" ", end = '')
    print(bcolors.ENDC+"\n")

    ##### DataGroup2 #####
    usr_in = input(f"optional DataGroup2 CSV ({len(DataGroup1)} elements required): ").replace(" ", "").replace("\t", "").replace("\n", "")
    if usr_in:
        val = re.sub(pattern, replace_range_by_ints, usr_in).split(",")
        for single in val:
            try:
                data_id = int(single) - 1
            except ValueError:
                print(bcolors.FAIL+'ERROR! Not a valid input!'+bcolors.ENDC)
                pause()
                exit(2)
            
            if data_id >= 0 and data_id <= len(data):
                if data_id not in DataGroup2:
                    DataGroup2.append(data_id)
                else:
                    print(bcolors.WARNING+"WARNING: '"+str(data_id+1)+"' is entered twice! Removing second instance!"+bcolors.ENDC)
            else:
                print(bcolors.WARNING+"WARNING: '"+str(data_id+1)+"' is out of range! Removing it!"+bcolors.ENDC)

        #print(DataGroup2)
        if len(DataGroup1) != len(DataGroup2):
            print(bcolors.FAIL+'ERROR! DataGroups not the same size!'+bcolors.ENDC)
            pause()
            exit(2)

        print(bcolors.OKBLUE+"DataGroup2: ", end = '')
        for key in DataGroup2:
            print(data[key]+" ", end = '')
        print(bcolors.ENDC+"\n")
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

suffixes = ['LFQ intensity ', 'Intensity ']
patterns = [re.compile(suffix) for suffix in suffixes]

def remove_obvious(s: str) -> str:
    for pattern in patterns:
        s = pattern.sub("", s)
    return s

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
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')

##############################################################################
#               Uporabniški vnosi
##############################################################################

#exp name
usr_name = slugify(input("Enter experiment name: "))
if not usr_name:
    usr_name = 'PeptideVisulizer'

now = datetime.now()
exp_name = now.strftime("%Y-%m-%d_%H-%M_") + usr_name
del now

print('Waiting for required files ... ')
# uporabnik izbere fajl peptides.txt in fasta fajl
peptides_file = askopenfilename(title = "Select peptides file",filetypes=[("peptides.txt", "peptides.txt")])
if not peptides_file:
    print(bcolors.FAIL+'ERROR! No peptides.txt file provided. Terminating...'+bcolors.ENDC)
    exit(2)
work_dir = os.path.dirname(peptides_file) + '/' + exp_name + '/'
fasta_file = askopenfilename(title = "Select fasta file",filetypes=[("FASTA file", ".fasta")])
if not fasta_file:
    print(bcolors.FAIL+'ERROR! No FASTA file provided. Terminating...'+bcolors.ENDC)
    exit(2)
# Y/n vprašanja
use_LFQ = query_yes_no('Do you want to use LFQ intensities?')
# peptide_threshold
try:
    peptide_threshold = float(input("Relative peptide detection threshold [default: {}%]: ".format(default_peptide_threshold)).replace(",", ".").replace(" ", "").replace("%", ""))
except ValueError:
    print(bcolors.WARNING+'WARNING! Relative peptide detection threshold could not be interpreted. Using default value '+str(default_peptide_threshold)+'%'+bcolors.ENDC)
    peptide_threshold = default_peptide_threshold
if peptide_threshold < 0 or peptide_threshold > 100:
    print(bcolors.WARNING+'WARNING! Relative peptide detection threshold makes no sense ('+str(peptide_threshold)+'%). Using default value '+str(default_peptide_threshold)+'%'+bcolors.ENDC)
    peptide_threshold = default_peptide_threshold

peptide_threshold *= 0.01

start_time = time.time()

# Replicates
try:
    replicates = int(input("Number of replicates used [default: 1]: ").replace(" ", ""))
except ValueError:
    print(bcolors.WARNING+'WARNING! Number of replicates could not be interpreted. Using default value 1'+bcolors.ENDC)
    replicates = 1
if replicates < 1 or replicates > replicates_MAX:
    print(bcolors.WARNING+'WARNING! Number of replicates makes no sense ('+str(replicates)+'). Using default value 1'+bcolors.ENDC)
    replicates = 1

##############################################################################
#               FASTA handling          ToDo: improve! + Uniprot download
##############################################################################

print('Reading FASTA ... ')
FASTA = {}
with open(fasta_file) as raw_fasta_file:
    for line in raw_fasta_file:
        #print(line)
        if line[0] == '>':
            sp=line[1:].split(' ')
            name = sp[0]
            FASTA[name] = ''
        else:
            FASTA[name] += line.strip()

print('Reading and filtering peptides.txt ... ')
with open(peptides_file) as raw_peptides_file:
    csv_reader = csv.reader(raw_peptides_file, delimiter="\t")
    line_count = 0
    for row in csv_reader:
        ###
        ##### v prvi vrstici ugotovim kateri stolpci me zanimajo + uporabnik pove kateri zanimajo njega #####
        ###
        if line_count == 0:
            sequence_index = row.index('Sequence')
            protein_index = row.index('Leading razor protein')
            
            start_position_index = row.index('Start position')
            length_index = row.index('Length')
            
            contaminant_index = row.index('Potential contaminant')
            reverse_index = row.index('Reverse')
            #LFQ intensity/Intensity?
            if use_LFQ:
                data_titles = [i for i in row if i.startswith('LFQ intensity ')]
                if not data_titles:
                    data_titles = [i for i in row if i.startswith('Intensity ')]
                    print(bcolors.WARNING+'WARNING: LFQ data not found ... falling back to intensities...'+bcolors.ENDC)
            else:
                data_titles = [i for i in row if i.startswith('Intensity ')]
            
            if not data_titles:
                print(bcolors.FAIL+'ERROR! No data detected!'+bcolors.ENDC)
                pause()
                exit()
            
            data_indexes = []
            for data_title in data_titles:
                data_indexes.append(row.index(data_title))
                #print(f'{data_title} z ID {row.index(data_title)} == {data_indexes[ data_titles.index(data_title) ]}')
            line_count += 1
            
            #odstranim nepotrebne dodatke imen
            if use_LFQ:
                data_titles = [w.replace('LFQ intensity ', '') for w in data_titles]
            else:
                data_titles = [w.replace('Intensity ', '') for w in data_titles]
            
            # v DB bo ključ protein in peptidi seznam
            DB = {}
            # tukaj bo ključ peptid in intenzitete seznam
            PEPTIDES_quantitative = {}
            # tukaj bo ključ peptid in bodo napisani peptidi od kje do kje
            PEPTIDES_qualitative = {}
        ###
        ##### Pridobivanje podatkov #####
        ###
        else:
            if row[reverse_index] != '+' and row[contaminant_index] != '+':
                tmp_protein = row[protein_index]
                tmp_peptide = row[sequence_index]
                tmp_start = int(row[start_position_index]) - 1
                tmp_length = int(row[length_index])
                PEPTIDES_qualitative[tmp_peptide] = [tmp_start, tmp_length]
                
                tmp_data = []
                for index in data_indexes:
                    tmp_data.append(float(row[index]))
            
                DB.setdefault(tmp_protein, []).append(tmp_peptide)
                PEPTIDES_quantitative[tmp_peptide] = tmp_data
            
                #just in case
                del tmp_protein
                del tmp_peptide
                del tmp_start
                del tmp_length
                del tmp_data
            line_count += 1

read_time = time.time() - start_time

###
##### Uporabnik nastavi skupine in način normalizacije
###

DataGroup1, DataGroup2 = userGrouping(data_titles)

if DataGroup2:
    group_normalization = query_yes_no('Do you want to normalize peptide intensities in each group?', 'no')
    print("") #tolko da je nova vrstica
else:
    group_normalization = False

continue_time = time.time()
N_max = len(DataGroup1)
##############################################################################
#               Generiranje rezultatov
##############################################################################

###
##### Normalizacija taka ko mora bit
# generiranje SUMMARY[peptide] = [SUM intensity, # peptides]
# generiranje SLICE_quantitative[protein][slice] = [sum, error(?)]
###

print("Statistically analyzing quantities, normalizing it and calculating error ... ")
#total slice intensity + normalize 
SLICE_quantitative = {}
for this_protein, these_peptides in DB.items():
    SLICE_quantitative[this_protein] = {}
    slice = 1
    for ID in DataGroup1:
        slice_intensity = 0
        for this_peptide in these_peptides:
            slice_intensity += PEPTIDES_quantitative[this_peptide][ID]
        SLICE_quantitative[this_protein][slice] = slice_intensity
        slice += 1
    if DataGroup2:
        slice = N_max+1
        for ID in DataGroup2:
            slice_intensity = 0
            for this_peptide in these_peptides:
                slice_intensity += PEPTIDES_quantitative[this_peptide][ID]
            SLICE_quantitative[this_protein][slice] = slice_intensity
            slice += 1

    slice_normalization_factor = max(SLICE_quantitative[this_protein].values())
    if slice_normalization_factor > 0:
        for n in range(len(SLICE_quantitative[this_protein].keys())):
            SLICE_quantitative[this_protein][n+1] = SLICE_quantitative[this_protein][n+1] / slice_normalization_factor
    del slice_normalization_factor

SUMMARY = {}
if group_normalization:
    print(bcolors.WARNING+"WARNING! Group normalization is used insted of overall normalization!"+bcolors.ENDC)
    
    for this_protein, these_peptides in DB.items():
        
        #DataGroup1
        tmp_list=[]
        peptides_count = 0
        for this_peptide in these_peptides:
            for ID in DataGroup1:
                tmp_list.append(PEPTIDES_quantitative[this_peptide][ID])
            peptides_count+=1
        #SUMMARY
        SUMMARY[this_protein] = [slugify(this_protein.replace("|", "_")),sum(tmp_list),peptides_count]
        
        normalization_factor = max(tmp_list)
        if normalization_factor > 0:
            for this_peptide in these_peptides:
                for ID in DataGroup1:
                    PEPTIDES_quantitative[this_peptide][ID] = PEPTIDES_quantitative[this_peptide][ID] / normalization_factor
        del normalization_factor
        del tmp_list
        #DataGroup2
        tmp_list=[]
        for this_peptide in these_peptides:
            for ID in DataGroup2:
                tmp_list.append(PEPTIDES_quantitative[this_peptide][ID])
        #SUMMARY
        SUMMARY[this_protein][1] = sum(tmp_list, SUMMARY[this_protein][1])
        
        normalization_factor = max(tmp_list)
        if normalization_factor > 0:
            for this_peptide in these_peptides:
                for ID in DataGroup2:
                    PEPTIDES_quantitative[this_peptide][ID] = PEPTIDES_quantitative[this_peptide][ID] / normalization_factor
        del normalization_factor
        del tmp_list
    
else:
    if DataGroup2:
        normalization_group = DataGroup1 + DataGroup2
    else:
        normalization_group = DataGroup1

    for this_protein, these_peptides in DB.items():
        tmp_list=[]
        peptides_count = 0
        for this_peptide in these_peptides:
            for ID in normalization_group:
                tmp_list.append(PEPTIDES_quantitative[this_peptide][ID])
            peptides_count+=1
        #SUMMARY
        SUMMARY[this_protein] = [slugify(this_protein.replace("|", "_")),sum(tmp_list),peptides_count]
        
        normalization_factor = max(tmp_list)
        if normalization_factor > 0:
            #print(tmp_list," - faktor: ",normalization_factor,"\n\n")
            for this_peptide in these_peptides:
                for ID in normalization_group:
                    PEPTIDES_quantitative[this_peptide][ID] = PEPTIDES_quantitative[this_peptide][ID] / normalization_factor
                #print(PEPTIDES_quantitative[this_peptide])
        del normalization_factor
        del tmp_list

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
    print(f'Generating PROTOMAP results in "{work_dir}"  ...')
else:
    print(f'Generating SEQUENCE ALIGN results in "{work_dir}" ... ')

if not os.path.isdir(work_dir):
    os.mkdir(work_dir)

with alive_bar(len(DB.keys()), dual_line=True, theme='classic') as bar: #loading bar
    for protein,peptide_list in DB.items():
        bar.text = f'Visulizing protein {protein}'
        #SUMMARY
        total_coverage = [False] * len(FASTA[protein])
        #PRIPRAVIM FAJLE
        result_file = work_dir+slugify(protein.replace("|", "_"))+'.txt' #TXT
        result_figure = work_dir+slugify(protein.replace("|", "_"))+'.svg' #SVG
        if os.path.exists(result_file): #TXT
            os.remove(result_file)
        if os.path.exists(result_figure): #SVG
            os.remove(result_figure)
        
        result = codecs.open(result_file, "w", "utf-8") #TXT
        result.write(u'\ufeff')
        if protein in FASTA and len(FASTA[protein]) > 0:
            result.write("SEQ:\t"+FASTA[protein]+"\n") #ToDo: nastavit ustrezno število tabov
        else:
            print(bcolors.FAIL+"ERROR! Given fasta file doesn't match fasta results from MaxQuant!\nERROR found at entry '",protein,"'!"+bcolors.ENDC)
            pause()
            exit()
    
        # IMAGE GENERATION START
        protein_length = len(FASTA[protein])

        labels = [] #ToDo: user defined labels or MWs
        for i in range(N_max):
            if DataGroup2:
                labels.append(data_titles[DataGroup1[i]]+'/'+data_titles[DataGroup2[i]])
            else:
                labels.append(data_titles[DataGroup1[i]])
    
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5), gridspec_kw={'width_ratios': [2, 1.2]})
        fig.suptitle(protein)

        ax1.set_xlim(1, protein_length)
        ax1.set_ylim(N_max+0.5, 0.5)
        ax1.set_xticks([1, protein_length])
        ax1.set_yticks(list(range(1, N_max+1)), labels)
        ax1.set_ylabel('Slice') #ob nastavitvi setMW spremeniti v MW [kDa]
        ax1.set_xlabel('Amino acid position', labelpad=0)
        
        secAx1 = ax1.twiny()
        secAx1.set_xlim(1, protein_length)
        secAx1.set_xticks([1, protein_length], ["N-term", "C-term"])
        
        ax2.set_xlim(0, 1)
        ax2.set_ylim(N_max+0.5, 0.5)
        ax2.set_xlabel('Slice intensity', labelpad=10)
        ax2.tick_params(bottom=False)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
            
        #nastavitev zamikov
        label_length_delta = len(max(labels, key=len)) - 5
        if label_length_delta > 0:
            y_offest = 0.1 + (label_length_delta * 0.009) # ToDo: preveriti, je če cifra ok
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
            if DataGroup2:
                experiment = data_titles[DataGroup1[i]]+'/'+data_titles[DataGroup2[i]]
            else:
                experiment = data_titles[DataGroup1[i]]
            
            coverage_list_Group1 = list(empty_char * len(FASTA[protein]))
            if DataGroup2:
                coverage_list_Group2 = list(empty_char * len(FASTA[protein]))
            
            for peptide in peptide_list:
                peptide_start = PEPTIDES_qualitative[peptide][0]
                peptide_length = PEPTIDES_qualitative[peptide][1]
                
                peptide_intensity_Group1 = PEPTIDES_quantitative[peptide][DataGroup1[i]]
                if DataGroup2:
                    peptide_intensity_Group2 = PEPTIDES_quantitative[peptide][DataGroup2[i]]
                alpha_intensity_Group1 = (1-alpha_boost)*peptide_intensity_Group1 + alpha_boost
                if DataGroup2:
                    alpha_intensity_Group2 = (1-alpha_boost)*peptide_intensity_Group2 + alpha_boost
                
                if DataGroup2:
                    #DataGroup1
                    if peptide_intensity_Group1 > peptide_threshold:
                        # GRAPH(SVG) OUTPUT MAIN
                        ax1.add_patch(Rectangle((peptide_start, slice-(bar_height/2)), peptide_length, bar_height/2,alpha=alpha_intensity_Group1 , linewidth=0, color='#283272'))
                        # TXT RESULTS
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
                            print(bcolors.WARNING+'WARNING: peptide "',peptide, '" does not fit onto protein "', protein, '"! Disregarding it in "',experiment,'" (', PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1], ' > ', len(coverage_list_Group1), ')'+bcolors.ENDC)
                    elif peptide_intensity_Group1 > 0 and NOTICE:
                        print('NOTICE: peptide',peptide, 'falls out of', round(peptide_threshold*100, 0) ,'% threshold (',round(peptide_intensity_Group1*100, 3),'%). Disregarding it in ', protein, data_titles[DataGroup1[i]])
                    #DataGroup2
                    if peptide_intensity_Group2 > peptide_threshold:
                        # GRAPH(SVG) OUTPUT MAIN
                        ax1.add_patch(Rectangle((peptide_start, slice), peptide_length, bar_height/2,alpha=alpha_intensity_Group2 , linewidth=0, color='#7f1710'))
                        # TXT RESULTS
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
                            print(bcolors.WARNING+'WARNING: peptide "',peptide, '" does not fit onto protein "', protein, '"! Disregarding it in "',experiment,'" (', PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1], ' > ', len(coverage_list_Group2), ')'+bcolors.ENDC)
                    elif peptide_intensity_Group2 > 0 and NOTICE:
                        print('NOTICE: peptide',peptide, 'falls out of', round(peptide_threshold*100, 0) ,'% threshold (',round(peptide_intensity_Group2*100, 3),'%). Disregarding it in ', protein, data_titles[DataGroup2[i]])
                else:
                    #DataGroup1
                    if peptide_intensity_Group1 > peptide_threshold:
                        # GRAPH(SVG) OUTPUT MAIN
                        ax1.add_patch(Rectangle((peptide_start, slice-(bar_height/2)), peptide_length, bar_height, alpha=alpha_intensity_Group1 , linewidth=0, color='#283272'))
                        # TXT RESULTS
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
                            print(bcolors.WARNING+'WARNING: peptide "',peptide, '" does not fit onto protein "', protein, '"! Disregarding it in "',experiment,'" (', PEPTIDES_qualitative[peptide][0] + PEPTIDES_qualitative[peptide][1], ' > ', len(coverage_list_Group1), ')'+bcolors.ENDC)
                    elif peptide_intensity_Group1 > 0 and NOTICE:
                        print('NOTICE: peptide "',peptide, '" falls out of', round(peptide_threshold*100, 0) ,'% threshold (',round(peptide_intensity_Group1*100, 3),'% ). Disregarding it in ', protein, data_titles[DataGroup1[i]])
            # dodam posamezni slice v txt in narišem skupno intenziteto slica + napako, če je podana
            coverage_Group1 = ''.join(coverage_list_Group1)
            result.write(data_titles[DataGroup1[i]]+"\t"+coverage_Group1+"\n")
            del coverage_list_Group1
            if DataGroup2:
                coverage_Group2 = ''.join(coverage_list_Group2)
                result.write(data_titles[DataGroup2[i]]+"\t"+coverage_Group2+"\n")
                del coverage_list_Group2
                # rišem dve intenziteti
                if replicates > 1:
                    this_err_1 = SLICE_quantitative_err[protein][slice]
                    this_err_2 = SLICE_quantitative_err[protein][N_max+slice]
                else:
                    this_err_1 = this_err_2 = 0
                ax2.barh(slice-(bar_width/4), SLICE_quantitative[protein][slice], xerr=this_err_1, height=(bar_width/2), color='#283272', alpha=0.5)
                ax2.barh(slice+(bar_width/4), SLICE_quantitative[protein][N_max+slice], xerr=this_err_2, height=(bar_width/2), color='#7f1710', alpha=0.5)

            else:
                #rišem eno intenziteto
                if replicates > 1:
                    this_err = SLICE_quantitative_err[protein][slice]
                else:
                    this_err = 0.0
                ax2.barh(slice, SLICE_quantitative[protein][slice], xerr=this_err, height=bar_width, color='#283272', alpha=0.5)
            
            slice+=1
        #SUMMARY
        SUMMARY[protein].append(round((sum(total_coverage) / len(total_coverage))*100, 1))
        #TXT KONEC
        result.close()
        #SVG KONEC
        plt.savefig(result_figure, format='svg', dpi=1200)
        plt.close()
        bar()

###
##### Generiranje HTML datoteke z rezultati
###
filenames = []
#for protein,summary_data in SUMMARY.items():
for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][3], reverse=True):
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
			.main_colum {width:60%;padding-left:2em;}
			.sub_colum {width:20%;padding-left:2em;}
			
		</style>
		<script>
			const PROTEINS = '''+json.dumps(filenames)+''';
			N = PROTEINS.length - 1;
			last_call = -1;
			
			function loadProtein(ID) {
				document.getElementById("content").innerHTML = '<img src="'''+exp_name+'''/'+PROTEINS[ID]+'.svg" alt="'+PROTEINS[ID]+'">';
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
					<tr style="text-align:left;">
						<th class="main_colum">Protein</th>
                        <th class="sub_colum">% coverage</th>
						<th class="sub_colum"># peptides</th>
						<th class="sub_colum">Gross intensity</th>
					</tr>
				</thead>
				<tbody>
'''
n = 0
#for protein,summary_data in SUMMARY.items():
for protein,summary_data in sorted(SUMMARY.items(), key=lambda e: e[1][3], reverse=True):
    HTML_CODE += '''				    <tr onclick="loadProtein('''+str(n)+''')">
						<td>'''+protein+'''</td>
                        <td>'''+str(summary_data[3])+'''</td>
						<td>'''+str(summary_data[2])+'''</td>
						<td>'''+str(summary_data[1])+'''</td>
                  </tr>
'''
    n+=1

HTML_CODE += '''
				</tbody>
			</table>
		</div>
		<div id="content"><br>'''+"Experiment \""+usr_name+"\"<br>LFQ intensities used? "+str(use_LFQ)+"<br>Relative threshold: "+str((peptide_threshold * 100))+"%<br>"+groups_summary+'''
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
print(f"\n{bcolors.OKGREEN}DONE!{bcolors.ENDC}\n{bcolors.OKBLUE}Read {line_count} peptides in {round(read_time,2)} seconds and processed results in {round(process_time/60,1)} minutes.{bcolors.ENDC}")

if DEBUG:
    if os.path.exists("debug.txt"):
        os.remove("debug.txt")
    file = open("debug.txt", "w")
    file.write("peptide_threshold: ")
    file.write(str(peptide_threshold*100))
    file.write("%\ndata_titles: ")
    file.write(json.dumps(data_titles))
    file.write("\ndata_indexes: ")
    file.write(json.dumps(data_indexes))
    file.write("\nDataGroup1: ")
    file.write(json.dumps(DataGroup1))
    file.write("\nDataGroup2: ")
    file.write(json.dumps(DataGroup2))
    file.write("\nSUMMARY: ")
    file.write(json.dumps(SUMMARY))
    file.write("\nDB: ")
    file.write(json.dumps(DB))
    file.write("\nSLICE_quantitative: ")
    file.write(json.dumps(SLICE_quantitative))
    file.write("\nPEPTIDES_qualitative: ")
    file.write(json.dumps(PEPTIDES_qualitative))
    file.write("\nPEPTIDES_quantitative: ")
    file.write(json.dumps(PEPTIDES_quantitative))
    file.write("FASTA: ")
    file.write(json.dumps(FASTA))
    file.close()
pause()