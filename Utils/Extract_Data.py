#!/usr/bin/env python3
import sys
import os
import pwd

from Scope.Classes_Data import *
from Scope.Thermal_Corrections import *
from Scope.Workflow import Branch
#from Scope import Constants

def extract_thermal_data(sys: object, branch_keyword: str, Trange: range=range(10,501,1), overwrite: bool=False, debug: int=0):

    ##############
    ### BRANCH ###
    ##############
    exists, this_branch = sys.find_branch(branch_keyword, debug=debug)
    print("Branch loaded with keyword", this_branch.keyword)
    #print("Nrecipes:", len(this_branch.recipes))
    if not exists: return False

    ## Checks that both spin states have an associated geometry minimum
    minima = True
    for recipe in this_branch.recipes:
        if not hasattr(recipe.subject,"min_coord"): minima = False; print(f"{recipe.subject.spin} of {sys.refcode} is not a minimum")
       
    ##################################
    ### RECIPE/MOLECULE PROPERTIES ###
    ##################################
    ## Only if both minima exist, the results are extracted. First for each individual molecule
    if not minima: 
        print(f"{sys.refcode} didn't reach both minima") 
    else: 
        print("Both molecules are minima")
        for recipe in this_branch.recipes:
            gmol = recipe.subject 
            print("    Doing recipe with:")
            print("    ",recipe.keyword, gmol.spin)
            ## Helec ##
            if overwrite or not "Helec" in recipe.results.keys():
                recipe.add_result(data("Helec",gmol.Helec,'au',"extract_thermal_data"), overwrite=overwrite)
                #recipe.add_result(data("Helec",gmol.Helec*Constants.har2kJmol,'kj',"extract_thermal_data"), overwrite=overwrite)
                print(gmol.spin, "Helec")
                recipe.results["Helec"].format()
            ## Selec ##
            if overwrite or not "Selec" in recipe.results.keys():
                recipe.add_result(get_Selec(gmol.spin, outunits='au'), overwrite=overwrite)
                recipe.results["Selec"].format()
            ## Hvib and Svib ##
            if overwrite or not "Hvib" in recipe.results.keys():
                Hvib_collection = collection("Hvib")
                for temp in Trange:
                    Hvib_collection.add_data(get_Hvib(gmol.freqs_cm, temp, freq_units='cm', outunits='au'))
                recipe.add_result(Hvib_collection, overwrite=overwrite)
                recipe.results["Hvib"].format()
            if overwrite or not "Svib" in recipe.results.keys():
                Svib_collection = collection("Svib")
                for temp in Trange:
                    Svib_collection.add_data(get_Svib(gmol.freqs_cm, temp, freq_units='cm', outunits='au'))
                recipe.add_result(Svib_collection, overwrite=overwrite)
                recipe.results["Svib"].format()
            ## Gtot ##
            if overwrite or not "Gtot" in recipe.results.keys():
                Gtot_collection = collection("Gtot")
                for temp in Trange:
                    # Retrieve data (not value)
                    Helec = recipe.results["Helec"]
                    Selec = recipe.results["Selec"]
                    Hvib = Hvib_collection.find_value_with_property("temp", temp)
                    Svib = Svib_collection.find_value_with_property("temp", temp)
                    #print(Helec.units)
                    #print(Selec.units)
                    #print(Hvib.units)
                    #print(Svib.units)
                    assert Helec.units == Selec.units == Hvib.units == Svib.units 
                    key = "Gtot"
                    value = get_Gibbs(Helec.value, Hvib.value, Selec.value, Svib.value, temp)
                    units = Helec.units 
                    function = "extract_thermal_data"
                    new_data = data(key, value, units, function)
                    new_data.add_property("temp", temp, overwrite=overwrite)
                    Gtot_collection.add_data(new_data)
                recipe.add_result(Gtot_collection, overwrite=overwrite)
                recipe.results["Gtot"].format()

    #########################
    ### BRANCH PROPERTIES ###
    #########################
    if minima:
        for recipe in this_branch.recipes:
            if   recipe.subject.spin == "HS": hs_rec = recipe 
            elif recipe.subject.spin == "LS": ls_rec = recipe 
            
        ## dHelec
        if overwrite or not "dHelec" in this_branch.results.keys():
            assert hs_rec.results["Helec"].units == ls_rec.results["Helec"].units 
            key = "dHelec" 
            value = hs_rec.results["Helec"].value - ls_rec.results["Helec"].value 
            units = hs_rec.results["Helec"].units
            function = "extract_thermal_data"
            this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

        ## dSelec
        if overwrite or not "dSelec" in this_branch.results.keys():
            assert hs_rec.results["Selec"].units == ls_rec.results["Selec"].units 
            key = "dSelec" 
            value = hs_rec.results["Selec"].value - ls_rec.results["Selec"].value 
            units = hs_rec.results["Selec"].units
            function = "extract_thermal_data"
            this_branch.add_result(data(key,value,units,function), overwrite=overwrite) 

        ## dHvib
        if overwrite or not "dHvib" in this_branch.results.keys():
            assert hs_rec.results["Hvib"].units == ls_rec.results["Hvib"].units 
            dHvib_collection = substract_collections("dHvib", hs_rec.results["Hvib"], ls_rec.results["Hvib"], prop="temp")
            this_branch.add_result(dHvib_collection, overwrite=overwrite)
            
        ## dSvib
        if overwrite or not "dSvib" in this_branch.results.keys():
            assert hs_rec.results["Svib"].units == ls_rec.results["Svib"].units 
            dSvib_collection = substract_collections("dSvib", hs_rec.results["Svib"], ls_rec.results["Svib"], prop="temp")
            this_branch.add_result(dSvib_collection, overwrite=overwrite)

        ## dGtot
        if overwrite or not "dGtot" in this_branch.results.keys():
            assert hs_rec.results["Gtot"].units == ls_rec.results["Gtot"].units 
            dGtot_collection = substract_collections("dGtot", hs_rec.results["Gtot"], ls_rec.results["Gtot"], prop="temp")
            this_branch.add_result(dGtot_collection, overwrite=overwrite)
            this_branch.results["dGtot"].format()

        ## T12
            print(f"Finding T12 between {Trange[0]}K and {Trange[-1]}K") 
            dGdata = dGtot_collection.get_values()
            print(f"dG between {dGdata[0]} and {dGdata[-1]}") 
            key = "T12"
            value = find_t12(Trange, dGdata)
            print(f"T12 value= {value}")
            #value = find_t12(Trange, dGtot_collection.datas)
            units = "K"
            function = "extract_thermal_data"
            this_branch.add_result(data(key, value, units, function), overwrite=overwrite)

        if debug > 0:
            print("Printing RESULTS for branch")
            this_branch.results["dHelec"].format()
            this_branch.results["dSelec"].format()
            this_branch.results["dHvib"].format()
            this_branch.results["dSvib"].format()
            this_branch.results["dGtot"].format()
            this_branch.results["T12"].format()

    return True
