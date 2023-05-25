import os
from ast import literal_eval

#######################
##### INPUT CLASS #####
#######################
class input_data(object):
    def __init__(self, f_name: str='', section: str='', debug=0):
        if f_name != '': self.read(f_name, section, debug=debug)
 
    def read(self, f_name: str, section: str, debug: int=0):
        dct = dict()
        with open(f_name, 'r') as f:
            canread = False
            for idx, line in enumerate(f.readlines()):
                if debug > 0: print("Doing line", idx, line) 

                if line[0] == '#' or line == '\n': continue
                line = line.strip()

                # Establishes the section of the input to read
                if line.startswith(section) and not canread: 
                    if debug > 0: print("Start of section", section, "found in line", idx) 
                    canread = True #; continue
                elif line.startswith('/') and canread:
                    canread = False #; continue
                    if debug > 0: print("end of section", section, "found in line", idx) 
                
                if canread:
                    data        = line.split('#')[0]
                    if len(data.split('=')) == 2: 
                        key, value  = data.split('=')
                        key         = key.strip()
                        value       = value.strip()
                        if debug > 0: print("key:", key)
                        if debug > 0: print("value:", value)
                        dct[key] = value
        self.set(dct)
        return self

    def set(self, dct):
        self.dct = dict()
        for k in dct.keys():
            attr        = self._set_attr(k, dct[k])
            self.dct[k] = attr
        return self

    def _set_attr(self, key, value):
        try:
            attr = literal_eval(value)
        except:
            attr = value
        setattr(self, key, attr)
        return attr

    def __repr__(self):
        string    = 'self.{:15}| {:20}| {:10}\n'
        to_print  = 'Formatted input interpretation: ( self -> Instance of class Input() )\n'
        to_print += '---------------------------------------------------\n'
        to_print += string.format('Key', 'Data Type', 'Value')
        to_print += '---------------------------------------------------\n'
        for key in self.dct.keys():
            val = self.dct[key]
            to_print += string.format(key, str(type(val)), str(val))
        return to_print

def set_resources_data(file_path, section="&resources", debug: int=0):
    resources = input_data(f_name=file_path, section=section, debug=debug)
    return resources

def set_options_data(file_path, section="&options", debug: int=0):
    options = input_data(f_name=file_path, section=section, debug=debug)
    return options

def set_job_data(file_path, section="&job_data", debug: int=0):
    job_data  = input_data(f_name=file_path, section=section, debug=debug)
    return job_data

def set_qc_data(file_path, section="&qc_data", debug: int=0):
    qc_data   = input_data(f_name=file_path, section=section, debug=debug)
    return qc_data

def interpret_software(name: str):
    name = name.lower()
    if name == 'g16' or name == 'gaussian': software = 'g16'
    elif name == 'qe' or name == 'quantum_espresso' or name == 'quantum_espresso' or name == 'espresso': software = 'qe'
    else:
        print("software",name," could not be interpreted")
        software = ''
    return software
