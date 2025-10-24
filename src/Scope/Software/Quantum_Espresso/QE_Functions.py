#from scope.classes_cell import Cell

################################################
## Functions Related to Quantum Espresso Only ##
################################################
def get_QE_data(state: object, debug: int=0):
    assert state.type == "state"
    if not hasattr(state,"moleclist"): state.get_moleclist()
    pairs = []
    for mol in state.moleclist:
        for at in mol.atoms:
            if tuple([at.label,at.spin]) not in pairs: pairs.append(tuple([at.label,at.spin]))
    if debug > 0: print(f"GET_QE_DATA. Pairs of label-spin found: {pairs}")
    return pairs

#####
def get_pp(elem: str, path: str):
    import json
    try: 
        with open(path+"config.json", 'r') as cfg:
            data = json.load(cfg)
            elemdata = data[elem]
            pp_name = elemdata["filename"]
            cutoff_wfc = elemdata["cutoff_wfc"]
            cutoff_rho = elemdata["cutoff_rho"]
    except Exception as exc: 
        raise ValueError(f"GET_PP: error reading JSON or parsing element {elem}. Printing Exception \n {exc}")
    return pp_name, cutoff_wfc, cutoff_rho
