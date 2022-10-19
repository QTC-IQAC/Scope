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

import Write_QE_Inputs
from Write_QE_Inputs import *

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

                                    gen_QE_iso_input(gmol,pathname,gmolname,"_scf",".inp",PPLIBRARY,'scf',isHubbard=True,isGrimme=True,U=2.27,cutoff=35,cubeside=70.0)
                                    gen_QE_subfile(pathname,gmolname,"_scf",".sub", "16", "iqtc09", "pw.x")

                                    gen_QE_iso_input(gmol,pathname,gmolname,"_relax",".inp",PPLIBRARY,'relax',isHubbard=True,isGrimme=True,U=2.27,cutoff=35,cubeside=70.0)
                                    gen_QE_subfile(pathname,gmolname,"_relax",".sub", "32", "iqtc09", "pw.x")
     
                                elif gmol.type == "Ligand":
                                    jobfolder=str(all_systems+'/'+core+'/'+crystal+'/LIGANDS/')

                                    gen_QE_iso_input(gmol,pathname,gmolname,"_scf",".inp",PPLIBRARY,'scf',isHubbard=True,isGrimme=True,U=2.27,cutoff=35,cubeside=70.0)
                                    gen_QE_subfile(pathname,gmolname,"_scf",".sub", "8", "iqtc09", "pw.x")

                                    gen_QE_iso_input(gmol,pathname,gmolname,"_relax",".inp",PPLIBRARY,'relax',isHubbard=True,isGrimme=True,U=2.27,cutoff=35,cubeside=70.0)
                                    gen_QE_subfile(pathname,gmolname,"_relax",".sub", "32", "iqtc09", "pw.x")
