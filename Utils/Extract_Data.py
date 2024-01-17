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

def extract_dH_solid(sys: object, branch_keyword: str, High_E_state: object, Low_E_state: object, overwrite: bool=False, global_env: object=None, debug: int=0):

    ## 0-Changes paths if necessary
    if global_env is not None:
        if global_env.check_paths(debug=1): 
            if debug > 1: print(f"EXECUTE_JOB, step 3b: global environment found with correct paths")
            try: 
                updated = sys.reset_paths(global_env, debug=0)
                if updated and debug > 1: print(f"EXECUTE_JOB, step 3b: system paths reset")
            except Exception as exc: 
                pass

    ### 1-Branch is loaded
    exists, this_branch = sys.find_branch(branch_keyword, debug=debug)
    if exists: print("Branch loaded with keyword", this_branch.keyword)
    if not exists: return False

    ### 2-Checks that energies have been parsed
    assert "energy" in High_E_state.results.keys()
    assert "energy" in Low_E_state.results.keys()

    ### 3-Verify properties of the state. Necessary for point 4
    if not hasattr(High_E_state,"fragmented"): High_E_state.check_fragmentation(reconstruct=True, debug=debug)
    if not hasattr(Low_E_state,"fragmented"): Low_E_state.check_fragmentation(reconstruct=True, debug=debug)
    assert not High_E_state.fragmented, f"Fragmented molecules in the geometry of High_E_state"
    assert not Low_E_state.fragmented, f"Fragmented molecules in the geometry of Low_E_state"
     
    ## 4-Get the number of complexes in the unit cells: this is why we need point 3
    ncomplex1 = 0
    for mol in High_E_state.moleclist:
        if mol.type == "Complex": ncomplex1 += 1
    ncomplex2 = 0
    for mol in Low_E_state.moleclist:
        if mol.type == "Complex": ncomplex2 += 1
        #if mol.type == "Complex" and hasattr(mol,"scope_guess_spin"): ncomplex2 += 1

    print("STEP4", ncomplex1, "complexes found in High_E_state")
    print("STEP4", ncomplex2, "complexes found in High_E_state")

    ## 5-Store Helec per molecule ##
    if overwrite or not "Helec" in High_E_state.results.keys():
        key = "Helec"
        units = High_E_state.results["energy"].units
        value = High_E_state.results["energy"].value / ncomplex1 
        High_E_state.add_result(data(key,value,str(units+"/molec"),"extract_dH_solid"), overwrite=overwrite)
    if overwrite or not "Helec" in Low_E_state.results.keys():
        key = "Helec"
        units = Low_E_state.results["energy"].units
        value = Low_E_state.results["energy"].value / ncomplex2 
        Low_E_state.add_result(data(key,value,str(units+"/molec"),"extract_dH_solid"), overwrite=overwrite)

    ## 6-Compute dHelec and store it in branch
    if overwrite or not "dHelec" in this_branch.results.keys():
        assert High_E_state.results["Helec"].units == Low_E_state.results["Helec"].units 
        key = "dHelec" 
        value = High_E_state.results["Helec"].value - Low_E_state.results["Helec"].value 
        units = High_E_state.results["Helec"].units
        function = "extract_dH_solid"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    if "dHelec" in this_branch.results.keys(): 
        this_branch.remove_output_lines()

    return True, this_branch.results["dHelec"]


def extract_T12(sys: object, branch_keyword: str, High_E_state: object, Low_E_state: object, Trange: range=range(10,501,1), flexible: bool=True, overwrite: bool=False, global_env: object=None, debug: int=0):

    ## 0-Changes paths if necessary
    if global_env is not None:
        if global_env.check_paths(debug=1): 
            if debug > 1: print(f"EXECUTE_JOB, step 3b: global environment found with correct paths")
            try: 
                updated = sys.reset_paths(global_env, debug=0)
                if updated and debug > 1: print(f"EXECUTE_JOB, step 3b: system paths reset")
            except Exception as exc: 
                pass
 
    if flexible:
        if not hasattr(High_E_state,"num_neq_freqs"): High_E_state.set_VNMs(High_E_state.VNMs)
        if not hasattr(Low_E_state,"num_neq_freqs"):  Low_E_state.set_VNMs(Low_E_state.VNMs)

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
 
    if not "Helec" in High_E_state.results: High_E_state.thermal_data()
    if not "Helec" in Low_E_state.results: Low_E_state.thermal_data()

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

    # Set branch_status_finished
    if "T12" in this_branch.results.keys(): 
        this_branch.set_status("finished", debug=debug)
        this_branch.remove_output_lines()

    return True, this_branch.results["T12"]
