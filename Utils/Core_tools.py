#!/usr/bin/env python3
import sys
import os
import pwd

def get_status(calcs_path: str, core: str, branch_keyword: str, debug: int=0):
    if calcs_path[-1] != '/' : calcs_path += '/'

    cancontinue = False
    if not os.path.isfile(calcs_path+core+'/'+core+".sys"):
        if debug > 1: print(f"System file of {core} not found in {calcs_path+core+'/'+core+'.sys'}")
    else:
        if os.path.isfile(calcs_path+core+'/'+"TERMINATED") and debug > 0:     print("Terminated", core)
        elif branch_keyword.lower() == "isolated" and os.path.isfile(calcs_path+core+'/'+"ISO_FINISHED") and debug > 0: print("Finished", core)
        elif branch_keyword.lower() == "solid" and os.path.isfile(calcs_path+core+'/'+"SOLID_FINISHED") and debug > 0: print("Finished", core)
        else: cancontinue = True
    return cancontinue

