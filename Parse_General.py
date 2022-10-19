import numpy as np

def search_string(string: str, lines: list, typ: str='all'):
    list_of_results = []
    for lx, line in enumerate(lines):
        if string in line:
            list_of_results.append(lx)  
    if len(list_of_results) > 0:
        if typ.lower() == 'all': return list(list_of_results), True
        elif typ.lower() == 'first': return list_of_results[0], True
        elif typ.lower() == 'last': return list_of_results[-1], True
        else: print("Inorrect variable 'typ' in search_string call")
    else: return [], False

#def search_string(string: str, lines: list, typ: str='all'):
#    list_of_results = []
#    for lx, line in enumerate(lines):
#        if string in line:
#            list_of_results.append(lx)  
#    if typ.lower() == 'all': return list(list_of_results)
#    elif typ.lower() == 'first': return list_of_results[0]
#    elif typ.lower() == 'last': return list_of_results[-1]
#    else: print("Correct variable 'typ' in search_string call")

def read_lines_file(filepath: str, flat: bool=False):
    info = open(filepath, 'r')
    lines = info.readlines()
    if flat: 
        for idx, l in enumerate(lines):
            lines[idx] = l.strip('\n')
    return np.array(lines)
