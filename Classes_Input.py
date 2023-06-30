import os
from ast import literal_eval

######################
def interpret_software(name: str):
    if name == 'g16' or name == 'gaussian': software = 'g16'
    elif name == 'qe' or name == 'quantum_espresso' or name == 'quantum_espresso' or name == 'espresso': software = 'qe'
    else: print("INTERPRET SOFTWARE: software",name," could not be interpreted. EXITING"); exit()
    return software

#######################
##### INPUT CLASS #####
#######################
class input_data(object):
    def __init__(self, f_name: str, section=None, debug=0):
        if f_name != '': self.read(f_name, section, debug=debug)
 
    def read(self, f_name: str, section=None, debug: int=0):
        dct = dict()
        with open(f_name, 'r') as f:
            canread = False
            for idx, line in enumerate(f.readlines()):
                if debug > 0: print("Doing line", idx, line) 

                if line[0] == '#' or line == '\n': continue
                line = line.strip()

                # Establishes the section of the input to read
                if section is not None:
                    if line.startswith(section) and not canread: 
                        if debug > 0: print("Start of section", section, "found in line", idx) 
                        canread = True #; continue
                    elif line.startswith('/') and canread:
                        canread = False #; continue
                        if debug > 0: print("end of section", section, "found in line", idx) 
                else: canread = True
                
                if canread:
                    data        = line.split('#')[0]
                    if len(data.split('=')) == 2: 
                        key, value  = data.split('=')
                        key         = key.strip()
                        value       = value.strip()
                        ### All capital letters are converted to lower letters
                        if type(key) == str:   key   = key.lower()   
                        if type(value) == str: 
                            if "/" not in value and key != 'branch': value  = value.lower() # Except paths
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
        try:      attr = literal_eval(value)
        except:   attr = value
        setattr(self, key, attr)
        return attr

    def _add_attr(self, dct, key, attr):
        self.dct[key] = attr
 
    def __repr__(self):
        string    = 'self.{:15}| {:20}| {:10}\n'
        to_print  = 'Formatted input interpretation: ( self -> Instance of class Input() )\n'
        to_print+= '---------------------------------------------------\n'
        to_print += string.format('Key', 'Data Type', 'Value')
        to_print += '---------------------------------------------------\n'
        for key in self.dct.keys():
            val = self.dct[key]
            to_print += string.format(key, str(type(val)), str(val))
        return to_print

#######################
def fill_environment_data(data: object, debug: int=0):
    ## Adds defaults to environment data
    return data

def fill_options_data(data: object, debug: int=0):
    ## Adds defaults to options data
    return data

def fill_job_data(data: object, debug: int=0):
    ## Adds defaults to job_data
    if not hasattr(data,"branch"):        print("WARNING: job_data is missing branch"); exit() 
    if not hasattr(data,"target"):        print("WARNING: job_data is missing target"); exit() 
    if not hasattr(data,"hierarchy"):     print("WARNING: job_data is missing hierarchy"); exit() 
    if not hasattr(data,"suffix"):        data._add_attr(suffix, str(data.hierarchy))
    if not hasattr(data,"keyword"):       data._add_attr(keyword, str(data.suffix))
    if not hasattr(data,"setup"):         data._add_attr(setup, "regular")
    if not hasattr(data,"requisites"):    data._add_attr(requisites, [])
    if not hasattr(data,"constrains"):    data._add_attr(constrains, ['self'])
    if not hasattr(data,"must_be_good"):  data._add_attr(must_be_good, False)
    return data

def fill_qc_data(data: object, debug: int=0):
    ## Adds defaults to qc_data
    if not hasattr(data,"software"):      print("WARNING: qc_data is missing software"); exit() 
    else:                                 data._add_attr(software, interpret_software(data.software))

    if data.software == "g16":
        if not hasattr(data,"functional"):    data._add_attr(functional, "B3LYP**")
        if not hasattr(data,"basis"):         data._add_attr(basis, "def2SVP")
        if not hasattr(data,"jobtype"):       data._add_attr(jobtype, "scf")
        if not hasattr(data,"loose_opt"):     data._add_attr(loose_opt, False)
        if not hasattr(data,"tight_opt"):     data._add_attr(tight_opt, False)
        if not hasattr(data,"is_grimme"):     data._add_attr(is_grimme, False)
        if not hasattr(data,"coord_tag"):     data._add_attr(coord_tag, "coord")

    elif data.software == "qe":
        if not hasattr(data,"coord_tag"):     data._add_attr(coord_tag, "coord")
        if not hasattr(data,"jobtype"):       data._add_attr(jobtype, "scf")
        if not hasattr(data,"functional"):    data._add_attr(functional, "pbe")
        if not hasattr(data,"is_hubbard"):    data._add_attr(is_hubbard, False)
        if not hasattr(data,"is_grimme"):     data._add_attr(is_grimme, False)
        if not hasattr(data,"uterm"): 
            if data.is_hubbard:               data._add_attr(uterm, float(2.27))
            else:                             data._add_attr(uterm, None)
        if not hasattr(data,"print_forces"):  data._add_attr(print_forces, False)
        if not hasattr(data,"cutoff"): 
            if data.jobtype == "scf":         data._add_attr(cutoff, int(25))
            elif data.jobtype == "relax":     data._add_attr(cutoff, int(25))
            elif data.jobtype == "vc-relax":  data._add_attr(cutoff, int(60))
        if not hasattr(data,"pressure"):      data._add_attr(pressure, int(0))
    
    return data

#######################

def set_environment_data(file_path, section="&environment", debug: int=0):
    environment = input_data(f_name=file_path, section=section, debug=debug)
    environment = fill_environment_data(environment)
    return environment

def set_options_data(file_path, section="&options", debug: int=0):
    options = input_data(f_name=file_path, section=section, debug=debug)
    options = fill_options_data(options)
    return options

def set_job_data(file_path, section="&job_data", debug: int=0):
    job_data = input_data(f_name=file_path, section=section, debug=debug)
    job_data = fill_job_data(job_data)
    return job_data

def set_qc_data(file_path, section="&qc_data", debug: int=0):
    qc_data = input_data(f_name=file_path, section=section, debug=debug)
    qc_data = fill_qc_data(qc_data)
    return qc_data

