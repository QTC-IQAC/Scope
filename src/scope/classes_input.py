import os
import sys
import numpy as np
from ast import literal_eval

######################
def interpret_software(name: str):
    if   name == 'g16' or name == 'gaussian': software = 'g16'
    elif name == 'qe' or name == 'quantum_espresso' or name == 'quantum_espresso' or name == 'espresso': software = 'qe'
    else: print("INTERPRET SOFTWARE: software",name," could not be interpreted. EXITING"); sys.exit()
    return software

#######################
##### INPUT CLASS #####
#######################
class Input_data(object):
    def __init__(self, content: str=None, section=None, isfile: bool=True, debug=0):
        
        if isfile: 
            content = os.path.abspath(content) 
            if not os.path.isfile(os.path.abspath(content)): raise FileNotFoundError(f"File not found: {content}")

        # Correct the section, in case the user forgot the &
        if section is not None:
            if section[0] != '&': section = f"&{section}"

        self.type     = "input_data"
        self.section  = section
        self.read(content, section, isfile=isfile, debug=debug)
 
    def read(self, source: str, section: str=None, isfile: bool=True, debug: int=0):
        dct = dict()

        ## Reads from file or content string
        if isfile: 
            from scope.parse_general import read_lines_file 
            lines    = read_lines_file(os.path.abspath(source)) 
        else:
            lines    = source.strip().split('\n')

        ## Reads Section Part
        if section is not None:      section_lines = read_section(lines, section=section, debug=0) 
        else:                        section_lines = lines
        if len(section_lines) == 0:  
            print(f"INPUT_DATA.READ: Section '{section}' is empty. Only defaults will be taken") 
            self.set({})
            return self
        elif debug > 0: print(f"INPUT_DATA.READ: {len(section_lines)} lines found in section {section}") 

        ## Creates Dictionary
        for idx, line in enumerate(section_lines):
            data        = line.split('#')[0]
            if len(data.split('=')) == 2: 
                key, value      = data.split('=')
                key             = key.strip()
                value           = value.strip()
                ### All capital letters are converted to lower letters
                if type(key)    == str:   key   = key.lower()   
                if type(value)  == str: 
                    if "/" not in value and key != 'branch' and value != 'True' and value != 'False': value  = value.lower() # Except paths
                if debug > 0: print("INPUT_DATA.READ: read entry -> key:", key, "value:", value)
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
        string    = 'self.{:20}| {:20}| {:10}\n'
        to_print  = 'Formatted input interpretation: ( self -> Instance of class Input() )\n'
        to_print+= '---------------------------------------------------\n'
        to_print += string.format('Key', 'Data Type', 'Value')
        to_print += '---------------------------------------------------\n'
        for key in self.dct.keys():
            val = self.dct[key]
            if key != 'section' and key != 'type':
                to_print += string.format(key, str(type(val)), str(val))
        return to_print

    def __eq__(self, other):
        if not isinstance(other, Input_data): return False
        else: 
            same = True
            for d in dir(self): 
                if '_' not in d and not isinstance(getattr(self, d), dict) and not callable(getattr(self, d)): 
                    at1 = getattr(self,d)
                    try:    
                        at2 = getattr(other,d)
                        if at1 != at2: 
                            same = False
                    except: return False
            for d in dir(other): 
                if '_' not in d and not isinstance(getattr(other, d), dict) and not callable(getattr(other, d)): 
                    at1 = getattr(other,d)
                    try:    
                        at2 = getattr(self,d)
                        if at1 != at2: 
                            same = False
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
    if not hasattr(data,"want_submit"):       data._add_attr("want_submit",True) 
    if not hasattr(data,"ignore_submitted"):  data._add_attr("ignore_submitted",False) 
    if not hasattr(data,"overwrite_inputs"):  data._add_attr("overwrite_inputs",True) 
    if not hasattr(data,"overwrite_outputs"): data._add_attr("overwrite_outputs",True) 
    return data

def fill_job_data(data: object, debug: int=0):
    ## Adds defaults to job_data
    if not hasattr(data,"branch"):        raise ValueError("WARNING: job_data is missing 'branch' input variable")
    if not hasattr(data,"hierarchy"):     raise ValueError("WARNING: job_data is missing 'hierarchy' input variable")  ## This shouldn't be mandatory, I need to fix it
    if not hasattr(data,"workflow"):      data._add_attr("workflow",str('all')),    
    if not hasattr(data,"suffix"):        data._add_attr("suffix", str(data.hierarchy))
    if not hasattr(data,"keyword"):       data._add_attr("keyword", str(data.suffix))
    if not hasattr(data,"istate"):        data._add_attr("istate", str("initial"))
    if not hasattr(data,"fstate"):        data._add_attr("fstate", str(data.suffix))
    if not hasattr(data,"job_setup"):     data._add_attr("job_setup", "regular")
    if not hasattr(data,"requisites"):    data._add_attr("requisites", [])
    if not hasattr(data,"constrains"):    data._add_attr("constrains", ['self'])
    if not hasattr(data,"must_be_good"):  data._add_attr("must_be_good", False)

    ## Adds defaults for rep_opt job setup type
    if data.job_setup == 'rep_opt' and not hasattr(data,"energy_thres"): data._add_attr("energy_thres",float(1e-5))
    if data.job_setup == 'rep_opt' and not hasattr(data,"max_steps"):    data._add_attr("max_steps",int(10))
    if data.job_setup == 'rep_opt' and not hasattr(data,"energiess"):    data._add_attr("energies",np.zeros((data.max_steps)))

    ## Modifies some attributes to avoid blank spaces and dashes, and to use lower letters
    data._mod_attr("keyword",str(data.keyword.lower().replace("-","_").replace(" ","_")))
    data._mod_attr("istate",str(data.istate.lower().replace("-","_").replace(" ","_")))
    data._mod_attr("fstate",str(data.fstate.lower().replace("-","_").replace(" ","_")))

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
    if not hasattr(data,"software"):      raise ValueError("ERROR: qc_data is missing 'software' input variable")
    else:                                 data._add_attr("software", interpret_software(data.software))

    if data.software == "g16":
        if not hasattr(data,"functional"):    data._add_attr("functional", "b3lyp**")
        if not hasattr(data,"basis"):         data._add_attr("basis", "sto-3g")
        if not hasattr(data,"jobtype"):       data._add_attr("jobtype", "scf")
        if not hasattr(data,"loose_opt"):     data._add_attr("loose_opt", False)
        if not hasattr(data,"tight_opt"):     data._add_attr("tight_opt", False)
        if not hasattr(data,"is_grimme"):     data._add_attr("is_grimme", False)
        if not hasattr(data,"grimme_type"):   data._add_attr("grimme_type", "d2")

    elif data.software == "qe":
        if not hasattr(data,"version"):       data._add_attr("version", float(7.0))
        if not hasattr(data,"pp_library"):    data._add_attr("pp_library", "vanderbilt")
        if not hasattr(data,"jobtype"):       data._add_attr("jobtype", "scf")
        if not hasattr(data,"functional"):    data._add_attr("functional", "pbe")
        if not hasattr(data,"is_hubbard"):    data._add_attr("is_hubbard", False)
        if not hasattr(data,"is_grimme"):     data._add_attr("is_grimme", False)
        if not hasattr(data,"grimme_type"):   data._add_attr("grimme_type", "d3bj")
        if not hasattr(data,"uterm"): 
            if data.is_hubbard:               data._add_attr("uterm", float(2.35)) ## Was 2.27 for some reason
            else:                             data._add_attr("uterm", None)
        if not hasattr(data,"print_forces"):  data._add_attr("print_forces", False)
        if not hasattr(data,"mix_beta"):      data._add_attr("mix_beta", float(0.7)) 
        if not hasattr(data,"elec_maxstep"):  data._add_attr("elec_maxstep", int(200)) 
        if not hasattr(data,"pressure"):      data._add_attr("pressure", int(0))
        if not hasattr(data,"forc_conv"):     data._add_attr("forc_conv", float(1e-5))
        if not hasattr(data,"elec_conv"):     data._add_attr("elec_conv", float(1e-5))

        available_jobtypes = ["opt", "relax", "freq", "scf", "vc-relax"]
        if data.jobtype not in available_jobtypes: raise ValueError(f"{data.jobtype} is not implemented")

    return data

#######################
def set_environment_data(content, section="&environment", isfile: bool=True, debug: int=0):
    environment = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    environment = fill_environment_data(environment)
    return environment

def set_options_data(content, section="&options", isfile: bool=True, debug: int=0):
    options = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    options = fill_options_data(options)
    return options

def set_job_data(content, section="&job_data", isfile: bool=True, debug: int=0):
    job_data = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    job_data = fill_job_data(job_data)
    return job_data

def set_qc_data(content, section="&qc_data", isfile: bool=True, debug: int=0):
    qc_data = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    qc_data = fill_qc_data(qc_data)
    return qc_data

#######################
def read_section(lines: list, section: str, debug: int=0):
    # Establishes the section of the input to read
    section_lines = []
    canread = False
    for idx, line in enumerate(lines):
        if line.startswith(section) and not canread: 
            if debug > 0: print("INPUT_DATA.READ_SECTION: Start of section", section, "found in line", idx) 
            canread = True #; continue
        elif line.startswith('/') and canread:
            canread = False #; continue
            if debug > 0: print("INPUT_DATA.READ_SECTION: end of section", section, "found in line", idx) 
        elif not line.startswith('/') and canread and len(line) > 0:
            if debug > 0: print(f"INPUT_DATA.READ_SECTION: reading line {line.strip()}")
            if line[0] == '#' or line == '\n': continue  ### Ignores lines commented with #
            line = line.strip()
            section_lines.append(line)
    return section_lines
