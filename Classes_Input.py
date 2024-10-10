import os
from ast import literal_eval
from Scope.Environment import set_cluster, set_user

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
        self.type = "input_data"
 
    def read(self, f_name: str, section=None, debug: int=0):
        dct = dict()
        with open(f_name, 'r') as f:
            canread = False
            for idx, line in enumerate(f.readlines()):

                if line[0] == '#' or line == '\n': continue
                line = line.strip()
                if debug > 0: print("INPUT DATA: Doing line", idx, line) 

                # Establishes the section of the input to read
                if section is not None:
                    if line.startswith(section) and not canread: 
                        if debug > 0: print("INPUT_DATA: Start of section", section, "found in line", idx) 
                        canread = True #; continue
                    elif line.startswith('/') and canread:
                        canread = False #; continue
                        if debug > 0: print("INPUT_DATA: end of section", section, "found in line", idx) 
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
                            if "/" not in value and key != 'branch' and value != 'True' and value != 'False': value  = value.lower() # Except paths
                        if debug > 0: print("key:", key, "value:", value)
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

    def _add_attr(self, key: str, value):
        try:      attr = literal_eval(value)
        except:   attr = value
        self.dct[key] = attr
        setattr(self, key, attr)

    def _mod_attr(self, key: str, value):
        try:      attr = literal_eval(value)
        except:   attr = value
        self.dct[key] = attr
        setattr(self, key, attr)
 
    def __add__(self, other):
        if not isinstance(other, type(self)): return self
        for d in dir(other):
            if '_' not in d and not callable(getattr(other,d)) and d != "dct":
                at1 = getattr(other,d)
                self._add_attr(d,at1)
        return self

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

    def __eq__(self, other):
        if not isinstance(other, input_data): return False
        else: 
            same = True
            for d in dir(self): 
                if '_' not in d and not callable(getattr(self,d)):
                    at1 = getattr(self,d)
                    try:    
                        at2 = getattr(other,d)
                        if at1 != at2: same = False
                    except: return False
            for d in dir(other): 
                if '_' not in d and not callable(getattr(other,d)): 
                    at1 = getattr(other,d)
                    try:    
                        at2 = getattr(self,d)
                        if at1 != at2: same = False
                    except: return False
        return same

#######################
def fill_environment_data(data: object, debug: int=0):
    ## Adds defaults to environment data
    if not hasattr(data,"requested_procs"):   data._add_attr("requested_procs",int(1)) 
    if not hasattr(data,"max_jobs"):          data._add_attr("max_jobs",int(100)) 
    if not hasattr(data,"max_procs"):         data._add_attr("max_procs",int(320)) 
    if not hasattr(data,"method"):            data._add_attr("method","weighted") 
    return data

def fill_options_data(data: object, debug: int=0):
    ## Adds defaults to options data
    return data

def fill_job_data(data: object, debug: int=0):
    ## Adds defaults to job_data
    if not hasattr(data,"branch"):        print("WARNING: job_data is missing branch"); exit() 
    if not hasattr(data,"target"):        print("WARNING: job_data is missing target"); exit() 
    if not hasattr(data,"hierarchy"):     print("WARNING: job_data is missing hierarchy"); exit() 
    if not hasattr(data,"suffix"):        data._add_attr("suffix", str(data.hierarchy))
    #if not hasattr(data,"keyword"):       data._add_attr("keyword", str("scf"))  ## In case I implement the available_keywords below
    if not hasattr(data,"keyword"):       data._add_attr("keyword", str(data.suffix))
    if not hasattr(data,"istate"):        data._add_attr("istate", str("initial"))
    if not hasattr(data,"fstate"):        data._add_attr("fstate", str(data.suffix))
    if not hasattr(data,"setup"):         data._add_attr("setup", "regular")
    if not hasattr(data,"requisites"):    data._add_attr("requisites", [])
    if not hasattr(data,"constrains"):    data._add_attr("constrains", ['self'])
    if not hasattr(data,"must_be_good"):  data._add_attr("must_be_good", False)
    #if not hasattr(data,"max_attempts") and data.must_be_good:  data._add_attr("max_attempts", int(5)) ## Not Implemented

    ## Modifies some attributes to avoid blank spaces and dashes, and to use lower letters
    data._mod_attr("keyword",str(data.keyword.lower().replace("-","_").replace(" ","_")))
    data._mod_attr("istate",str(data.istate.lower().replace("-","_").replace(" ","_")))
    data._mod_attr("fstate",str(data.fstate.lower().replace("-","_").replace(" ","_")))

   # available_keywords = ["opt", "relax", "freq", "scf", "vc-relax"]
   # if data.keyword not in available_keywords:
   #     print("----------------------------------")
   #     print("WARNING: job_keyword not available")
   #     print("----------------------------------"); exit() 

    ## Modifies the list of requisites and constrains to avoid blanks and dashes
    tmp_req = []
    for r in data.requisites:
        tmp_req.append(r.replace("-","_").replace(" ","_"))
    data._mod_attr("requisites",tmp_req)
    tmp_con = []
    for c in data.constrains:
        tmp_con.append(c.replace("-","_").replace(" ","_"))
    data._mod_attr("constrains",tmp_con)

    return data

def fill_qc_data(data: object, debug: int=0):
    ## Adds defaults to qc_data
    if not hasattr(data,"software"):      print("WARNING: qc_data is missing software"); exit() 
    else:                                 data._add_attr("software", interpret_software(data.software))

    if data.software == "g16":
        if not hasattr(data,"functional"):    data._add_attr("functional", "B3LYP**")
        if not hasattr(data,"basis"):         data._add_attr("basis", "def2SVP")
        if not hasattr(data,"jobtype"):       data._add_attr("jobtype", "scf")
        if not hasattr(data,"loose_opt"):     data._add_attr("loose_opt", False)
        if not hasattr(data,"tight_opt"):     data._add_attr("tight_opt", False)
        if not hasattr(data,"is_grimme"):     data._add_attr("is_grimme", False)

    elif data.software == "qe":
        if not hasattr(data,"jobtype"):       data._add_attr("jobtype", "scf")
        if not hasattr(data,"functional"):    data._add_attr("functional", "pbe")
        if not hasattr(data,"is_hubbard"):    data._add_attr("is_hubbard", False)
        if not hasattr(data,"is_grimme"):     data._add_attr("is_grimme", False)
        if not hasattr(data,"uterm"): 
            if data.is_hubbard:               data._add_attr("uterm", float(2.35)) ## Was 2.27 for some reason
            else:                             data._add_attr("uterm", None)
        if not hasattr(data,"print_forces"):  data._add_attr("print_forces", False)
        if not hasattr(data,"cutoff"): 
            if data.jobtype == "scf":         data._add_attr("cutoff", int(25))
            elif data.jobtype == "relax":     data._add_attr("cutoff", int(25))
            elif data.jobtype == "vc-relax":  data._add_attr("cutoff", int(60))
        if not hasattr(data,"mix_beta"):      data._add_attr("mix_beta", float(0.7)) 
        if not hasattr(data,"elec_maxstep"):  data._add_attr("elec_maxstep", int(200)) 
        if not hasattr(data,"pressure"):      data._add_attr("pressure", int(0))
    
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

