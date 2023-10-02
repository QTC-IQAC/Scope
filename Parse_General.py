import numpy as np

def search_string(string: str, lines: list, typ: str='all', lowlim: int=None, uplim: int=None, debug: int=0):
    if type(lines) != list: lines = list([lines])
    if debug > 0: print(f"SEARCH STRING: searching {string} in file")
    list_of_results = []
    for lx, line in enumerate(lines):
        if string in line:
            if debug > 0: print(f"SEARCH STRING: found string in line {lx}")
            if lowlim and uplim:
                if lx >= lowlim and lx <= uplim: list_of_results.append(lx)  
            elif lowlim and not uplim:
                if lx >= lowlim: list_of_results.append(lx)  
            elif uplim and not lowlim:
                if lx <= uplim: list_of_results.append(lx)  
            elif not lowlim and not uplim:
                list_of_results.append(lx)  
    if len(list_of_results) > 0:
        if typ.lower() == 'all': return list(list_of_results), True
        elif typ.lower() == 'first': return list_of_results[0], True
        elif typ.lower() == 'last': return list_of_results[-1], True
        else: print("Incorrect variable 'typ' in search_string call")
    else: return [], False


def read_lines_file(filepath: str, flat: bool=False):
    info = open(filepath, 'r')
    lines = info.readlines()
    if flat: 
        for idx, l in enumerate(lines):
            lines[idx] = l.strip('\n')
    return np.array(lines)

def slurm_time_to_seconds(sl_time: str):
    if sl_time == "infinite":                 ## Infinite 
        days = 100
        time = days * 86400 
    elif ':' in sl_time and '-' in sl_time:   ## Standard Format
        blocks = sl_time.replace('-',' ').replace(':',' ').split()
        if len(blocks) == 4: 
            days     = int(blocks[0])
            hours    = int(blocks[1])
            minutes  = int(blocks[2])
            seconds  = int(blocks[3])
            time = days * 86400 + hours * 3600 + minutes * 60 + seconds
        elif len(blocks) == 3: 
            days     = int(blocks[0])
            hours    = int(blocks[1])
            minutes  = int(blocks[2])
            seconds  = 0 
            time = days * 86400 + hours * 3600 + minutes * 60 + seconds
        else: print(f"slurm_time_to_seconds: slurm time could not be parsed: {sl_time}")
    elif ':' in sl_time and '-' not in sl_time:   ## No days
        blocks = sl_time.replace(':',' ').split()
        if len(blocks) == 3: 
            days     = 0 
            hours    = int(blocks[0])
            minutes  = int(blocks[1])
            seconds  = int(blocks[2])
            time = days * 86400 + hours * 3600 + minutes * 60 + seconds
        else: print(f"slurm_time_to_seconds: slurm time could not be parsed: {sl_time}")
    else:
        print(f"slurm_time_to_seconds: slurm time could not be parsed: {sl_time}")
        time = 0
    return int(time)
