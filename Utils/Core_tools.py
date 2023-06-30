#!/usr/bin/env python3
import sys
import os
import pwd

def get_status(sys_path: str, core: str, branch_keyword, debug: int=0):
    ## Keep in mind that here, the system is not yet loaded
    if sys_path[-1] != '/' : sys_path += '/'

    status = "active"
    if not os.path.isfile(sys_path+core+'/'+core+".sys"):
        if debug > 1: print(f"System file of {core} not found in {sys_path+core+'/'+core+'.sys'}")
        return None

    else:
        if os.path.isfile(sys_path+core+'/'+"TERMINATED") and debug > 0:   
            status = "terminated"
        else:
            if type(branch_keyword) == str:
                if branch_keyword.lower() == "isolated" and os.path.isfile(sys_path+core+'/'+"ISO_FINISHED"):
                    status = "finished"
                elif branch_keyword.lower() == "solid" and os.path.isfile(sys_path+core+'/'+"SOLID_FINISHED"): 
                    status = "finished"
            elif type(branch_keyword) == list:
                sk = []
                for bk in branch_keyword:
                    if bk.lower() == "isolated":
                        if os.path.isfile(sys_path+core+'/'+"ISO_FINISHED"): sk.append("finished")
                        else:                                                sk.append("active")
                    elif bk.lower() == "solid":
                        if os.path.isfile(sys_path+core+'/'+"SOLID_FINISHED"): sk.append("finished")
                        else:                                                  sk.append("active")
                if "active" in sk: status = "active"
                else:              status = "finished"
            else: print("GET_STATUS: I could not understand branch_keyword:", branch_keyword)
    if debug > 0: print(core, "status:", status)
    return str(status)

