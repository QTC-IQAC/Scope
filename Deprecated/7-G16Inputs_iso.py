import pickle
import sys
import os
import shutil
import csv

from cell2mol.tmcharge_common import Cell, atom, molecule, ligand, metal
from cell2mol.readwrite import writexyz, print_molecule
from cell2mol.tmcharge_common import labels2formula
from cell2mol.elementdata import ElementData
from cell2mol.readwrite import readcif, search_string_in_file

utilspath = '/home/g4vela/SCOPE/Database_SCO/Scripts'
sys.path.append(utilspath)

import Write_G16_Inputs
from Write_G16_Inputs import *

all_systems = os.getcwd()

PPLIBRARY = '/home/g4vela/SCOPE/Database_SCO/Scripts/renamed_PP_Library/'

list_of_csv_entries = []
for core in sorted(os.listdir(all_systems)):
    #print("doing core:", core)
    if os.path.isdir(all_systems+'/'+core):
        for crystal in sorted(os.listdir(all_systems+'/'+core)):
            #print("doing crystal:", crystal)
            if os.path.isdir(all_systems+'/'+core+'/'+crystal):
                for fil in os.listdir(all_systems+'/'+core+'/'+crystal):
                    #print("doing fil:", fil)
                    if fil.endswith(".gmol") and core in fil:
 
                        pathfile=all_systems+'/'+core+'/'+crystal+'/'+fil
                        #pathname=pathfile.split(".")[0]
                        pathname=all_systems+'/'+core+'/'+crystal+'/'
                        gmolname=fil.split(".")[0]
                        with open(pathfile, "rb") as fil:
                            gmol = pickle.load(fil)
                        
                            if hasattr(gmol,"type"):
                                if gmol.type == "Complex": 
                                    jobfolder=str(all_systems+'/'+core+'/'+crystal+'/TMCs/')

                                    gen_G16_iso_input(gmol,pathname,gmolname,"_opt_HS",".com",jobtype='opt&freq',functional='B3LYP*',basis='def2SVP', spin='HS', isGrimme=True, nproc=16)
                                    gen_G16_subfile(pathname,gmolname, suffix="_opt_HS", extension=".sub", cores="16", queue="iqtc09")

                                    gen_G16_iso_input(gmol,pathname,gmolname,"_opt_LS",".com",jobtype='opt&freq',functional='B3LYP*',basis='def2SVP', spin='LS', isGrimme=True, nproc=16)
                                    gen_G16_subfile(pathname,gmolname, suffix="_opt_LS", extension=".sub", cores="16", queue="iqtc09")

                                elif gmol.type == "Ligand":
                                    jobfolder=str(all_systems+'/'+core+'/'+crystal+'/LIGANDS/')

                                    gen_G16_iso_input(gmol,pathname,gmolname,"_opt",".com",jobtype='opt&freq',functional='B3LYP*',basis='def2SVP',isGrimme=True, nproc=16)
                                    gen_G16_subfile(pathname,gmolname, suffix="_opt", extension=".sub", cores="8", queue="iqtc09")
