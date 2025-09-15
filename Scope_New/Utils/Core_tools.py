#!/usr/bin/env python3
import sys
import os
import pwd

def get_status(sys_path: str, core: str, branch_keyword, debug: int=0):
    if sys_path[-1] != '/' : sys_path += '/'
    ##########################################################
    ## Keep in mind that here, the system is not yet loaded ##
    ## The status is set based, only, on the files present  ##
    ##########################################################
    ## We do not want to load the system to save time       ##
    ##########################################################

    status = "active"

    ## If system file is not present
    if not os.path.isfile(sys_path+core+'/'+core+".sys"):
        if debug > 1: print(f"System file of {core} not found in {sys_path+core+'/'+core+'.sys'}")
        return "absent"

    ## Otherwise..
    if os.path.isfile(sys_path+core+'/'+"TERMINATED"):  return "terminated"

    ## If one branch keyword is given
    if type(branch_keyword) == str:
        if branch_keyword.lower() == "isolated" and os.path.isfile(sys_path+core+'/'+"ISO_FINISHED"):  return "finished"
        elif branch_keyword.lower() == "solid" and os.path.isfile(sys_path+core+'/'+"SOLID_FINISHED"): return "finished" 

    ## If more than one branch keyword is given
    elif type(branch_keyword) == list:
        sk = []
        for bk in branch_keyword:
            if bk.lower() == "isolated":
                if os.path.isfile(sys_path+core+'/'+"ISO_FINISHED"): sk.append("finished")
                else:                                                sk.append("active")
            elif bk.lower() == "solid":
                if os.path.isfile(sys_path+core+'/'+"SOLID_FINISHED"): sk.append("finished")
                else:                                                  sk.append("active")
        if "active" in sk: return "active"
        else:              return "finished"

    else: print("GET_STATUS: I could not understand branch_keyword:", branch_keyword)
    return "active"

