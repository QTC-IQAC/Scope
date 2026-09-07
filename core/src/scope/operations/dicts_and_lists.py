import numpy as np

####################
### Dictionaries ###
####################
def same_dictionaries(dic1: dict, dic2: dict):
    from collections import Counter
    """
    Checks if two specific types of dictionaries, the signatures obtained in scope.operations.graphs.get_signatures() are equivalent:
    """
    # 1) Compare that they have the same layers
    if not list(dic1.keys()) == list(dic2.keys()): return False
    # 2) Compares the counter
    for layer in list(dic1.keys()):
        if not Counter(dic1[layer].values()) == Counter(dic2[layer].values()): return False
    return True

#############
### Lists ###
#############
def extract_from_list(entrylist: list, old_array: list, dimension: int=2, debug: int=0) -> list:
    if debug >= 1: print(f"EXTRACT_FROM_LIST. received: {entrylist=}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. received: {old_array=}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. maximum value received in entrylist: {np.max(entrylist)+1}")
    if debug >= 1: print(f"EXTRACT_FROM_LIST. length of old_array: {len(old_array)}")
    assert len(old_array) >= np.max(entrylist)+1
    length = len(entrylist)
    if dimension == 2:
        new_array = np.empty((length, length), dtype=object)
        for idx, row in enumerate(entrylist):
            for jdx, col in enumerate(entrylist):
                new_array[idx, jdx] = old_array[row][col]
    elif dimension == 1:
        new_array = np.empty((length), dtype=object)
        for idx, val in enumerate(entrylist):
            new_array[idx] = old_array[val]
    return list(new_array)

def where_in_array(array,condition) -> list:
    results = []
    for idx, a in enumerate(array):
        if a == condition: results.append(idx)
    return results

def mergelists(list1, list2, prop1, prop2) -> list:
    #print("Received", list1, list2)
    nitems=len(list1)+len(list2)
    mergedlist = []
    for idx in range(0,nitems):
        for jdx, at1 in enumerate(list1):
            if (idx == at1):
                mergedlist.append(prop1[jdx])
        for jdx, at2 in enumerate(list2):
            if (idx == at2):
                mergedlist.append(prop2[jdx])
    return mergedlist

def range2list(rang: range) -> list:
    lst = []
    for i in rang:
        lst.append(i)
    return lst