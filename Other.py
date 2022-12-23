import sys

def where_in_array(array,condition) -> list:
    results = []
    for idx, a in enumerate(array):
        if a == condition: results.append(idx)
    return results

