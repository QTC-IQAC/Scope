#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Data import *
from Scope.Thermal_Corrections import *
from Scope import Constants

######
def get_T12(branch: object, High_E_state: object, Low_E_state: object, Trange: range=range(10,501,1), ignore_not_minima: bool=False, flexible: bool=True, overwrite: bool=False, debug: int=0):

    if not ignore_not_minima:
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
 
    assert "Helec" in High_E_state.results 
    assert "Helec" in Low_E_state.results 

    ## dHelec
    if overwrite or not "dHelec" in branch.results.keys():
        assert High_E_state.results["Helec"].units == Low_E_state.results["Helec"].units 
        key = "dHelec" 
        value = High_E_state.results["Helec"].value - Low_E_state.results["Helec"].value 
        units = High_E_state.results["Helec"].units
        function = "Scope.Spin_Crossover.SCO_functions.get_T12"
        branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dSelec
    if overwrite or not "dSelec" in branch.results.keys():
        assert High_E_state.results["Selec"].units == Low_E_state.results["Selec"].units 
        key = "dSelec" 
        value = High_E_state.results["Selec"].value - Low_E_state.results["Selec"].value 
        units = High_E_state.results["Selec"].units
        function = "Scope.Spin_Crossover.SCO_functions.get_T12"
        branch.add_result(data(key,value,units,function), overwrite=overwrite) 

    ## dHvib
    if overwrite or not "dHvib" in branch.results.keys():
        assert High_E_state.results["Hvib"].units == Low_E_state.results["Hvib"].units 
        dHvib_collection = substract_collections("dHvib", High_E_state.results["Hvib"], Low_E_state.results["Hvib"], prop="temperature")
        branch.add_result(dHvib_collection, overwrite=overwrite)
        
    ## dSvib
    if overwrite or not "dSvib" in branch.results.keys():
        assert High_E_state.results["Svib"].units == Low_E_state.results["Svib"].units 
        dSvib_collection = substract_collections("dSvib", High_E_state.results["Svib"], Low_E_state.results["Svib"], prop="temperature")
        branch.add_result(dSvib_collection, overwrite=overwrite)

    ## dGtot
    if overwrite or not "dGtot" in branch.results.keys():
        assert High_E_state.results["Gtot"].units == Low_E_state.results["Gtot"].units 
        dGtot_collection = substract_collections("dGtot", High_E_state.results["Gtot"], Low_E_state.results["Gtot"], prop="temperature")
        branch.add_result(dGtot_collection, overwrite=overwrite)

    ## T12
    if overwrite or not "T12" in branch.results.keys():
        dGdata = dGtot_collection.get_values()
        key = "T12"
        value = find_t12(Trange, dGdata)
        units = "K"
        function = "Scope.Spin_Crossover.SCO_functions.get_T12"
        branch.add_result(data(key, value, units, function), overwrite=overwrite)

    if debug > 0:
        print("Printing RESULTS for branch")
        branch.results["dHelec"].print_in_units('kj')
        print(branch.results["dSelec"])
        print(branch.results["dHvib"])
        print(branch.results["dSvib"])
        print(branch.results["dGtot"])
        print(branch.results["T12"])

    ## Set branch_status_finished
    #if "T12" in branch.results.keys(): 
    #    branch.set_status("finished", debug=debug)
    #    branch.remove_output_lines()

    return True, branch.results["T12"]
