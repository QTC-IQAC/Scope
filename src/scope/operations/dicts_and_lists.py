import numpy as np

####################
### Dictionaries ###
####################


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