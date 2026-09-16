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
def extract_from_list(entrylist: list, old_array: list, dimension: int=2) -> list:
    """
    Extract selected entries from an atom-indexed sequence or square matrix.

    Parameters
    ----------
    entrylist:  Indices to retain from the array.
    old_array:  With ``dimension=1``, a sequence whose entries may themselves be vectors or objects, such as labels, coordinates, or Atom objects.
                With ``dimension=2``, a square matrix indexed by atom along both axes, such as an adjacency matrix.
    dimension:  1 selects entries along the first axis. 2 selects the same indices along both axes (x,y) of a square matrix.
    """
    assert dimension in (1, 2), f"Unsupported dimension: {dimension}"
    assert all(0 <= index < len(old_array) for index in entrylist), \
        "entrylist contains indices outside old_array"

    # Case of dimension == 2
    if dimension == 2:
        assert all(len(row) == len(old_array) for row in old_array), \
            "dimension=2 requires a square matrix"
        return [[old_array[row][col] for col in entrylist] for row in entrylist]

    # Case of dimension == 1
    return [old_array[index] for index in entrylist]

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