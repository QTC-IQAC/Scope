#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Data import *
from Scope.Thermal_Corrections import *
from Scope.Workflow import Branch
from Scope import Constants

##############################################################################
### Here is where the all protocols to extract system properties, (often) involving more than one system (eg. HS and LS), are collected
##############################################################################

def extract_T12(sys: object, branch_keyword: str, High_E_state: object, Low_E_state: object, Trange: range=range(10,501,1), flexible: bool=True, overwrite: bool=False, global_env: object=None, debug: int=0):

    ## Uses data stored in state-class objects. More information can be found in Classes_State

    ## 0-Changes paths if necessary
    if global_env is not None:
        if global_env.check_paths(debug=1): 
            if debug > 1: print(f"EXECUTE_JOB, step 3b: global environment found with correct paths")
            try: 
                updated = sys.reset_paths(global_env, debug=0)
                if updated and debug > 1: print(f"EXECUTE_JOB, step 3b: system paths reset")
            except Exception as exc: 
                pass
 
    ### Not sure why this was necessary. Removed
    #if flexible:
    #    if not hasattr(High_E_state,"num_neg_freqs"): High_E_state.set_VNMs(High_E_state.VNMs)
    #    if not hasattr(Low_E_state,"num_neg_freqs"):  Low_E_state.set_VNMs(Low_E_state.VNMs)

    ##############
    ### BRANCH ###
    ##############
    exists, this_branch = sys.find_branch(branch_keyword, debug=debug)
    if exists and debug > 0: print("Branch loaded with keyword", this_branch.keyword)
    if not exists: return False, None

    ## Checks that both states exist and are minima
    if not hasattr(High_E_state,"isminimum") or not hasattr(Low_E_state,"isminimum"): 
        if debug > 0: print(f"{High_E_state.name} and/or {Low_E_state.name} have not been evaluated as minima")
        return False, None

    ## If not, it evaluates if the states could be considered minima for the sake of T12 evaluation
    if not High_E_state.isminimum and (not High_E_state.almost_minimum or not flexible):
        if debug > 0: print(f"{High_E_state.name} cannot be used for Thermochemistry")
        return False, None
        
    if not Low_E_state.isminimum and (not Low_E_state.almost_minimum or not flexible):
        if debug > 0: print(f"{High_E_state.name} cannot be used for Thermochemistry")
        return False, None
 
    if not "Helec" in High_E_state.results: High_E_state.get_thermal_data()
    if not "Helec" in Low_E_state.results: Low_E_state.get_thermal_data()

    ## dHelec
    if overwrite or not "dHelec" in this_branch.results.keys():
        assert High_E_state.results["Helec"].units == Low_E_state.results["Helec"].units 
        key = "dHelec" 
        value = High_E_state.results["Helec"].value - Low_E_state.results["Helec"].value 
        units = High_E_state.results["Helec"].units
        function = "extract_T12"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dSelec
    if overwrite or not "dSelec" in this_branch.results.keys():
        assert High_E_state.results["Selec"].units == Low_E_state.results["Selec"].units 
        key = "dSelec" 
        value = High_E_state.results["Selec"].value - Low_E_state.results["Selec"].value 
        units = High_E_state.results["Selec"].units
        function = "extract_T12"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dHvib
    if overwrite or not "dHvib" in this_branch.results.keys():
        assert High_E_state.results["Hvib"].units == Low_E_state.results["Hvib"].units 
        dHvib_collection = substract_collections("dHvib", High_E_state.results["Hvib"], Low_E_state.results["Hvib"], prop="temp")
        this_branch.add_result(dHvib_collection, overwrite=overwrite)
        
    ## dSvib
    if overwrite or not "dSvib" in this_branch.results.keys():
        assert High_E_state.results["Svib"].units == Low_E_state.results["Svib"].units 
        dSvib_collection = substract_collections("dSvib", High_E_state.results["Svib"], Low_E_state.results["Svib"], prop="temp")
        this_branch.add_result(dSvib_collection, overwrite=overwrite)

    ## dGtot
    if overwrite or not "dGtot" in this_branch.results.keys():
        assert High_E_state.results["Gtot"].units == Low_E_state.results["Gtot"].units 
        dGtot_collection = substract_collections("dGtot", High_E_state.results["Gtot"], Low_E_state.results["Gtot"], prop="temp")
        this_branch.add_result(dGtot_collection, overwrite=overwrite)

    ## T12
    if overwrite or not "T12" in this_branch.results.keys():
        dGdata = dGtot_collection.get_values()
        key = "T12"
        value = find_t12(Trange, dGdata)
        units = "K"
        function = "extract_T12"
        this_branch.add_result(data(key, value, units, function), overwrite=overwrite)

    if debug > 0:
        print("Printing RESULTS for branch")
        this_branch.results["dHelec"].print_in_units('kj')
        print(this_branch.results["dSelec"])
        print(this_branch.results["dHvib"])
        print(this_branch.results["dSvib"])
        print(this_branch.results["dGtot"])
        print(this_branch.results["T12"])

    #if debug > 1:
    #    for t in Trange: 

    # Set branch_status_finished
    if "T12" in this_branch.results.keys(): 
        this_branch.set_status("finished", debug=debug)
        this_branch.remove_output_lines()

    return True, this_branch.results["T12"]
