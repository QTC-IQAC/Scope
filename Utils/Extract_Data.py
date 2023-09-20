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

def extract_dH_solid(sys: object, branch_keyword: str, state1: object, state2: object, overwrite: bool=False, debug: int=0):

    ### 1-Branch is loaded
    exists, this_branch = sys.find_branch(branch_keyword, debug=debug)
    if exists: print("Branch loaded with keyword", this_branch.keyword)
    if not exists: return False

    #### 2-Checks that both states exist and are minima
    #if not hasattr(state1,"isminimum") or not hasattr(state2,"isminimum"): 
    #    if debug > 0: print(f"{state1.name} and/or {state2.name} are not evaluated as minima")
    #    return False, None
    #if not state1.isminimum or not state2.isminimum: 
    #    if debug > 0: print(f"{state1.name} and/or {state2.name} are not minima")
    #    return False, None

    ### 2-Checks that energies have been parsed
    assert "energy" in state1.results.keys()
    assert "energy" in state2.results.keys()

    ### 3-Verify properties of the state. Necessary for point 4
    if not hasattr(state1,"fragmented"): state1.check_fragmentation(reconstruct=True, debug=debug)
    if not hasattr(state2,"fragmented"): state2.check_fragmentation(reconstruct=True, debug=debug)
    assert not state1.fragmented, f"Fragmented molecules in the geometry of state1"
    assert not state2.fragmented, f"Fragmented molecules in the geometry of state2"
     
    ## 4-Get the number of complexes in the unit cells: this is why we need point 3
    ncomplex1 = 0
    for mol in state1.moleclist:
        if mol.type == "Complex": ncomplex1 += 1
    ncomplex2 = 0
    for mol in state2.moleclist:
        if mol.type == "Complex": ncomplex2 += 1
        #if mol.type == "Complex" and hasattr(mol,"scope_guess_spin"): ncomplex2 += 1

    print("STEP4", ncomplex1, "complexes found in state1")
    print("STEP4", ncomplex2, "complexes found in state1")

    ## 5-Store Helec per molecule ##
    if overwrite or not "Helec" in state1.results.keys():
        key = "Helec"
        units = state1.results["energy"].units
        value = state1.results["energy"].value / ncomplex1 
        state1.add_result(data(key,value,str(units+"/molec"),"extract_dH_solid"), overwrite=overwrite)
    if overwrite or not "Helec" in state2.results.keys():
        key = "Helec"
        units = state2.results["energy"].units
        value = state2.results["energy"].value / ncomplex2 
        state2.add_result(data(key,value,str(units+"/molec"),"extract_dH_solid"), overwrite=overwrite)

    ## 6-Compute dHelec and store it in branch
    if overwrite or not "dHelec" in this_branch.results.keys():
        assert state1.results["Helec"].units == state2.results["Helec"].units 
        key = "dHelec" 
        value = state1.results["Helec"].value - state2.results["Helec"].value 
        units = state1.results["Helec"].units
        function = "extract_dH_solid"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    if "dHelec" in this_branch.results.keys(): 
        this_branch.remove_output_lines()

    return True, this_branch.results["dHelec"]


def extract_T12(sys: object, branch_keyword: str, state1: object, state2: object, Trange: range=range(10,501,1), overwrite: bool=False, debug: int=0):

    ##############
    ### BRANCH ###
    ##############
    exists, this_branch = sys.find_branch(branch_keyword, debug=debug)
    if exists and debug > 0: print("Branch loaded with keyword", this_branch.keyword)
    if not exists: return False, None

    ## Checks that both states exist and are minima
    if not hasattr(state1,"isminimum") or not hasattr(state2,"isminimum"): return False, None
    if not state1.isminimum or not state2.isminimum: 
        if debug > 0: print(f"{state1.name} and/or {state2.name} are not minima")
        return False, None

    if not "Helec" in state1.results: state1.thermal_data()
    if not "Helec" in state2.results: state2.thermal_data()

    ## dHelec
    if overwrite or not "dHelec" in this_branch.results.keys():
        assert state1.results["Helec"].units == state2.results["Helec"].units 
        key = "dHelec" 
        value = state1.results["Helec"].value - state2.results["Helec"].value 
        units = state1.results["Helec"].units
        function = "extract_T12"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dSelec
    if overwrite or not "dSelec" in this_branch.results.keys():
        assert state1.results["Selec"].units == state2.results["Selec"].units 
        key = "dSelec" 
        value = state1.results["Selec"].value - state2.results["Selec"].value 
        units = state1.results["Selec"].units
        function = "extract_T12"
        this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dHvib
    if overwrite or not "dHvib" in this_branch.results.keys():
        assert state1.results["Hvib"].units == state2.results["Hvib"].units 
        dHvib_collection = substract_collections("dHvib", state1.results["Hvib"], state2.results["Hvib"], prop="temp")
        this_branch.add_result(dHvib_collection, overwrite=overwrite)
        
    ## dSvib
    if overwrite or not "dSvib" in this_branch.results.keys():
        assert state1.results["Svib"].units == state2.results["Svib"].units 
        dSvib_collection = substract_collections("dSvib", state1.results["Svib"], state2.results["Svib"], prop="temp")
        this_branch.add_result(dSvib_collection, overwrite=overwrite)

    ## dGtot
    if overwrite or not "dGtot" in this_branch.results.keys():
        assert state1.results["Gtot"].units == state2.results["Gtot"].units 
        dGtot_collection = substract_collections("dGtot", state1.results["Gtot"], state2.results["Gtot"], prop="temp")
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
