import os
import sys
import numpy as np
from ast import literal_eval
from fnmatch import fnmatchcase

######################
def interpret_software(name: str):
    if   name == 'g16' or name == 'gaussian': software = 'g16'
    elif name == 'qe' or name == 'quantum_espresso' or name == 'quantum_espresso' or name == 'espresso': software = 'qe'
    else: 
        raise ValueError(f"INTERPRET SOFTWARE: software {name} not recognized. Available options are 'g16' or 'qe'. EXITING")
    return software

#######################
##### INPUT CLASS #####
#######################
class Input_data(object):
    """
    Represents parsed key-value input data from a file or text block.

    Attributes:
        object_type (str):              Object category (`"input_data"`).
        section (str):                  Parsed section name.
        dct (dict):                     Parsed values as a dictionary.

    Methods:
        read():                         Parses raw content into key-value pairs.
        set():                          Stores a dictionary as attributes.
        _set_attr():                    Converts and assigns a parsed value.
        _add_attr():                    Adds a new parsed attribute.
        _mod_attr():                    Modifies an existing parsed attribute.
    """
    def __init__(self, content: str=None, section=None, isfile: bool=True, debug=0):
        
        if isfile: 
            if not os.path.isfile(content): 
                if os.path.isfile(os.path.abspath(content)): 
                    content = os.path.abspath(content) 
                else: raise FileNotFoundError(f"INPUT_DATA.INIT: Content File not found in {content} or {os.path.abspath(content)}")

        # Correct the section, in case the user forgot the &
        if section is not None:
            if section[0] != '&': section = f"&{section}"

        self.object_type = "input_data"
        self.section  = section
        self.read(content, section, isfile=isfile, debug=debug)
 
    def read(self, source: str, section: str=None, isfile: bool=True, debug: int=0):
        dct = dict()

        ## Reads from file or content string
        if isfile: 
            from scope.read_write import read_lines_file
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
        """Adds missing public settings from other, retaining existing values.

        Mutates and returns self; both + and += use this behavior. Existing
        None and False values are retained. Setting names may contain
        underscores; private attributes, methods and dct are excluded.
        """
        if not isinstance(other, type(self)): return self
        for d in dir(other):
            if not d.startswith('_') and d != 'dct' and not callable(getattr(other,d)) and not hasattr(self,d):
                self._add_attr(d,getattr(other,d))
        return self

    def __repr__(self):
        string    = 'self.{:20}| {:20}| {:10}\n'
        to_print  = 'Formatted input interpretation: ( self -> Instance of class Input() )\n'
        to_print+= '---------------------------------------------------\n'
        to_print += string.format('Key', 'Data Type', 'Value')
        to_print += '---------------------------------------------------\n'
        for key in self.dct.keys():
            val = self.dct[key]
            if key != 'section' and key != 'object_type':
                to_print += string.format(key, str(type(val)), str(val))
        return to_print

    def __eq__(self, other):
        """Compares public non-callable attributes, excluding duplicate dct storage.

        Public names may contain underscores and values may be dictionaries.
        NumPy array settings use exact element/shape equality. Section metadata
        remains part of the comparison.
        """
        if not isinstance(other, Input_data): return False
        else: 
            same = True
            for d in dir(self): 
                if not d.startswith('_') and d != 'dct' and not callable(getattr(self, d)):
                    at1 = getattr(self,d)
                    try:    
                        at2 = getattr(other,d)
                        if isinstance(at1, np.ndarray) or isinstance(at2, np.ndarray):
                            if not np.array_equal(at1, at2): same = False
                        elif at1 != at2:
                            same = False
                    except: return False
            for d in dir(other): 
                if not d.startswith('_') and d != 'dct' and not callable(getattr(other, d)):
                    at1 = getattr(other,d)
                    try:    
                        at2 = getattr(self,d)
                        if isinstance(at1, np.ndarray) or isinstance(at2, np.ndarray):
                            if not np.array_equal(at1, at2): same = False
                        elif at1 != at2:
                            same = False
                    except: return False
        return same

#######################
def check_environment_data(data: object, debug: int=0) -> None:
    """Validates supplied environment settings without changing data.

    Accepted values:
        - time_limit: positive integer in minutes or None ('none' in parsed input).

    Missing optional settings use defaults when filled. Raises ValueError for
    invalid data; returns None otherwise.
    """
    time_limit = getattr(data,'time_limit',None)
    if isinstance(time_limit,str) and time_limit == 'none': return
    if time_limit is not None and (isinstance(time_limit,bool) or not isinstance(time_limit,int) or time_limit <= 0):
        raise ValueError("CHECK_ENVIRONMENT_DATA: time_limit must be a positive integer in minutes, or None")

def check_options_data(data: object, debug: int=0) -> None:
    """Validates supplied submission/overwrite flags without changing data.

    Accepted values:
        - want_submit: True or False.
        - ignore_submitted: True or False.
        - overwrite_inputs: True or False.
        - overwrite_outputs: True or False.

    Missing flags use defaults when filled. Raises ValueError for invalid
    data; returns None otherwise.
    """
    for setting in ['want_submit','ignore_submitted','overwrite_inputs','overwrite_outputs']:
        if hasattr(data,setting) and not isinstance(getattr(data,setting),bool):
            raise ValueError(f"CHECK_OPTIONS_DATA: {setting} must be True or False")

def check_job_data(data: object, debug: int=0) -> None:
    """Validates required Job fields and supplied options without changing data.

    Accepted values:
        - job: nonempty string.
        - istate: nonempty string.
        - fstate: nonempty string; distinct from istate after name cleanup for 'findiff'.
        - recovery: nonempty string; strategy availability is checked by the software's recovery function.
        - save_wf: False, None, True, 'job' or 'workflow'; string aliases 'none', 'false' and 'true' are also accepted, ignoring case and surrounding spaces.
        - job_setup: nonempty string, e.g. 'regular', 'rep_opt' or 'findiff'.
        - max_steps: positive integer for 'rep_opt'.
        - findiff_displacement: finite positive number for 'findiff', in findiff_units.
        - findiff_units: 'bohr' or 'angstrom' for 'findiff'.

    Requires branch, workflow, job and hierarchy. Missing optional fields use
    defaults when filled. Raises ValueError for invalid data; returns None
    otherwise.
    """
    # Checks branch, workflow, job and hierarchy
    for setting in ['branch','workflow','job','hierarchy']:
        if not hasattr(data,setting): 
            raise ValueError(f"CHECK_JOB_DATA: job_data is missing {setting!r} input variable")

    # Checks job name, and istate/fstate name
    for setting in ['job','istate','fstate']:
        if hasattr(data,setting) and (not isinstance(getattr(data,setting),str) or not getattr(data,setting).strip()):
            raise ValueError(f"CHECK_JOB_DATA: {setting} must be a nonempty name")

    # Checks the recovery strategy name; the software checks its availability.
    recovery = getattr(data,'recovery','standard')
    if not isinstance(recovery,str) or not recovery.strip():
        raise ValueError("CHECK_JOB_DATA: recovery must be a nonempty strategy name")

    # Checks the wavefunction storage mode and its accepted aliases.
    save_wf = getattr(data,'save_wf',False)
    if isinstance(save_wf,str):
        if save_wf.strip().lower() not in ['none','false','true','job','workflow']:
            raise ValueError("CHECK_JOB_DATA: save_wf must be False, None, 'job' or 'workflow'")
    elif save_wf is not None and not isinstance(save_wf,bool):
        raise ValueError("CHECK_JOB_DATA: save_wf must be False, None, 'job' or 'workflow'")

    # Checks the setup name before evaluating its specific settings.
    job_setup = getattr(data,'job_setup','regular')
    if not isinstance(job_setup,str) or not job_setup.strip():
        raise ValueError("CHECK_JOB_DATA: job_setup must be a nonempty string")

    if job_setup == 'rep_opt':
        # Checks the number of repetitive optimization steps.
        max_steps = getattr(data,'max_steps',5)
        if isinstance(max_steps,bool) or not isinstance(max_steps,(int,np.integer)) or max_steps <= 0:
            raise ValueError("CHECK_JOB_DATA: max_steps must be a positive integer")
    elif job_setup == 'findiff':
        # Checks the finite-difference displacement length and units.
        displacement = getattr(data,'findiff_displacement',0.01)
        if isinstance(displacement,bool) or not isinstance(displacement,(int,float,np.integer,np.floating)) or not np.isfinite(displacement) or displacement <= 0:
            raise ValueError("CHECK_JOB_DATA: findiff_displacement must be finite and positive")
        if getattr(data,'findiff_units','bohr') not in ['bohr','angstrom']:
            raise ValueError("CHECK_JOB_DATA: findiff_units must be 'bohr' or 'angstrom'")

        # Checks that initial/final State names remain distinct after cleanup.
        istate = getattr(data,'istate','initial').strip().lower().replace('-','_').replace(' ','_')
        fstate = getattr(data,'fstate',data.job).strip().lower().replace('-','_').replace(' ','_')
        if fstate == istate:
            raise ValueError("CHECK_JOB_DATA: Finite differences require distinct reference and final States")

def check_qc_data(data: object, debug: int=0) -> None:
    """Validates the QC software and computation type without changing data.

    Accepted values:
        - software: 'g16' (alias 'gaussian') or 'qe' (aliases 'quantum_espresso' and 'espresso').
        - comp_type (G16): 'scf', 'opt', 'freq', 'ts', 'td' or 'tda'.
        - comp_type (QE): 'scf', 'opt', 'relax', 'vc-relax' or 'freq'.

    Requires software; missing comp_type means 'scf'. Raises ValueError for
    invalid data; returns None otherwise.
    """
    if not hasattr(data,'software'): raise ValueError("CHECK_QC_DATA: qc_data is missing 'software' input variable")
    software  = interpret_software(data.software)
    comp_type = getattr(data,'comp_type','scf')
    if software == 'g16': available_types = ['scf','opt','freq','ts','td','tda']
    elif software == 'qe': available_types = ['scf','opt','relax','vc-relax','freq']
    if comp_type not in available_types:
        raise ValueError(f"CHECK_QC_DATA: computation type {comp_type!r} is not implemented for {software}")

#######################
def fill_environment_data(data: object, debug: int=0):
    """Adds resource defaults and normalizes a parsed 'none' time limit.

    Expects data accepted by check_environment_data(); returns the same object.
    time_limit is in integer minutes; None uses the selected queue's limit.
    """
    ## Adds defaults to environment data
    if not hasattr(data,"requested_procs"):   data._add_attr("requested_procs",int(1)) 
    if not hasattr(data,"max_jobs"):          data._add_attr("max_jobs",int(100)) 
    if not hasattr(data,"max_procs"):         data._add_attr("max_procs",int(320)) 
    if not hasattr(data,"method"):            data._add_attr("method","weighted") 
    if not hasattr(data,"time_limit"):        data._add_attr("time_limit",None)
    
    if isinstance(data.time_limit,str) and data.time_limit == 'none': data._mod_attr('time_limit',None)
    return data

def fill_options_data(data: object, debug: int=0):
    """Adds submission/overwrite defaults and returns the same object.

    Expects data accepted by check_options_data().
    """
    ## Adds defaults to options data
    if not hasattr(data,"want_submit"):       data._add_attr("want_submit",True) 
    if not hasattr(data,"ignore_submitted"):  data._add_attr("ignore_submitted",False) 
    if not hasattr(data,"overwrite_inputs"):  data._add_attr("overwrite_inputs",True) 
    if not hasattr(data,"overwrite_outputs"): data._add_attr("overwrite_outputs",True) 
    return data

def fill_job_data(data: object, debug: int=0):
    """Adds Job defaults and normalizes names and accepted option aliases.

    Expects data accepted by check_job_data(); returns the same object.
    recovery defaults to 'standard'. save_wf=True becomes 'job', and None
    becomes False. Defaults for rep_opt/findiff depend on job_setup.
    """
    ## Adds defaults to job_data
    if not hasattr(data,"istate"):        data._add_attr("istate", str("initial"))
    if not hasattr(data,"fstate"):        data._add_attr("fstate", str(data.job))
    if not hasattr(data,"job_setup"):     data._add_attr("job_setup", "regular")
    if not hasattr(data,"requisites"):    data._add_attr("requisites", [])
    if not hasattr(data,"constrains"):    data._add_attr("constrains", ['self'])
    if not hasattr(data,"must_be_good"):  data._add_attr("must_be_good", True)
    if not hasattr(data,"save_wf"):       data._add_attr("save_wf", False)
    if not hasattr(data,"recovery"):      data._add_attr("recovery", "standard")

    data._mod_attr('recovery',data.recovery.strip().lower())

    # None also disables storage of the wavefunction. Input_data reads an unquoted None as the string 'none'.
    save_wf = data.save_wf
    if save_wf is None or save_wf is False: save_wf = False
    elif save_wf is True:                   save_wf = 'job'
    elif isinstance(save_wf, str):
        save_wf = save_wf.strip().lower()
        if save_wf in ['none', 'false']: save_wf = False
        elif save_wf == 'true':          save_wf = 'job'
    data._mod_attr('save_wf', save_wf)

    ## Adds defaults for rep_opt job setup type
    if data.job_setup == 'rep_opt' and not hasattr(data,"energy_thres"): data._add_attr("energy_thres",float(1e-5))
    if data.job_setup == 'rep_opt' and not hasattr(data,"max_steps"):    data._add_attr("max_steps",int(5))
    if data.job_setup == 'rep_opt' and not hasattr(data,"energies"):     data._add_attr("energies",np.zeros((data.max_steps)))

    # Finite-difference settings; geometry coordinates remain in Angstrom.
    if data.job_setup == 'findiff':
        if not hasattr(data,'findiff_displacement'): data._add_attr('findiff_displacement', 0.01)
        if not hasattr(data,'findiff_units'):        data._add_attr('findiff_units', 'bohr')

    ## Modifies some attributes to avoid blank spaces and dashes, and to use lower letters
    data._mod_attr("job",str(data.job.strip().lower().replace("-","_").replace(" ","_")))
    data._mod_attr("istate",str(data.istate.strip().lower().replace("-","_").replace(" ","_")))
    data._mod_attr("fstate",str(data.fstate.strip().lower().replace("-","_").replace(" ","_")))

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
    """Adds software/type-specific QC defaults and normalizes software aliases.

    Expects data accepted by check_qc_data(); returns the same object.
    """
    ## Adds defaults to qc_data
    data._mod_attr("software", interpret_software(data.software))

    if data.software == "g16":
        # Adds the default computation type.
        if not hasattr(data,"comp_type"):      data._add_attr("comp_type", "scf")

        # Adds defaults
        if not hasattr(data,"functional"):     data._add_attr("functional", "pbe")
        if not hasattr(data,"basis"):          data._add_attr("basis", "sto-3g")
        if not hasattr(data,"loose_opt"):      data._add_attr("loose_opt", False)
        if not hasattr(data,"tight_opt"):      data._add_attr("tight_opt", False)
        if not hasattr(data,"ultrafine_grid"): data._add_attr("ultrafine_grid", False)
        if not hasattr(data,"print_forces"):   data._add_attr("print_forces", False)
        if not hasattr(data,"unrestricted"):   data._add_attr("unrestricted", False)

        # Grimme Dispersion Corrections
        if not hasattr(data,"is_grimme"):      
            if not hasattr(data,"grimme_type"): data._add_attr("is_grimme", False)
            else:                               data._add_attr("is_grimme", True)
        if data.is_grimme: 
            if not hasattr(data,"grimme_type"): data._add_attr("grimme_type", "d3bj")

        if data.comp_type == "td" or data.comp_type == "tda":
            if not hasattr(data,"td_type"):    data._add_attr("td_type", "singlets")  ## Will be ignored if specie.spin > 0
            if not hasattr(data,"td_nstates"): data._add_attr("td_nstates", int(10))

        elif data.comp_type == "ts" or data.comp_type == "opt": 
            ## Default is recalcFC=30
            if not hasattr(data,"fctype"): data._add_attr("fctype","recalcfc")
            if data.fctype == "recalcfc" and not hasattr(data,"recalc_steps"): data._add_attr("recalc_steps", int(30))

    elif data.software == "qe":
        # Adds the default computation type.
        if not hasattr(data,"comp_type"):     data._add_attr("comp_type", "scf")

        # Adds defaults
        if not hasattr(data,"version"):       data._add_attr("version", float(7.0))
        if not hasattr(data,"pp_library"):    data._add_attr("pp_library", "vanderbilt")
        if not hasattr(data,"functional"):    data._add_attr("functional", "pbe")
        if not hasattr(data,"is_hubbard"):    data._add_attr("is_hubbard", False)
        if not hasattr(data,"is_grimme"):     data._add_attr("is_grimme", False)
        if not hasattr(data,"grimme_type"):   data._add_attr("grimme_type", "d3bj")
        if not hasattr(data,"uterm"): 
            if data.is_hubbard:               data._add_attr("uterm", float(2.35)) ## Was 2.27 for some reason
            else:                             data._add_attr("uterm", None)
        if not hasattr(data,"print_forces"):  data._add_attr("print_forces", False)
        if not hasattr(data,"mix_beta"):      data._add_attr("mix_beta", float(0.7)) 
        if not hasattr(data,"mix_mode"):      data._add_attr("mix_mode", "plain")
        if not hasattr(data,"mix_ndim"):      data._add_attr("mix_ndim", int(8))
        if not hasattr(data,"elec_maxstep"):  data._add_attr("elec_maxstep", int(200)) 
        if not hasattr(data,"opt_maxstep"):   data._add_attr("opt_maxstep", int(100))
        if not hasattr(data,"pressure"):      data._add_attr("pressure", int(0))
        if not hasattr(data,"forc_conv"):     data._add_attr("forc_conv", float(1e-4))
        if not hasattr(data,"elec_conv"):     data._add_attr("elec_conv", float(1e-8))
        if not hasattr(data,"trust_rad_ini"): data._add_attr("trust_rad_ini", float(0.1))
        if not hasattr(data,"trust_rad_min"): data._add_attr("trust_rad_min", float(1e-6))
        if not hasattr(data,"upscale"):       data._add_attr("upscale", int(100))
        if not hasattr(data,"bfgs_ndim"):     data._add_attr("bfgs_ndim", int(1))
        if not hasattr(data,"ion_dynamics"):  data._add_attr("ion_dynamics", "bfgs")
        if not hasattr(data,"cell_dynamics"): data._add_attr("cell_dynamics", "bfgs")

    return data

def harmonize_input_data(environment: object, options: object, job_data: object, qc_data: object, debug: int=0):
    if hasattr(job_data, "job_type"): job_data._mod_attr("job_type", qc_data.comp_type)
    else:                             job_data._add_attr("job_type", qc_data.comp_type)
    return environment, options, job_data, qc_data

#######################
def set_environment_data(content, section="&environment", isfile: bool=True, debug: int=0):
    """Parses and validates environment settings, then adds defaults."""
    environment = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    check_environment_data(environment, debug=max(debug-1,0))
    environment = fill_environment_data(environment)
    return environment

def set_options_data(content, section="&options", isfile: bool=True, debug: int=0):
    """Parses and validates submission options, then adds defaults."""
    options     = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    check_options_data(options, debug=max(debug-1,0))
    options     = fill_options_data(options)
    return options

def set_job_data(content, section="&job_data", isfile: bool=True, debug: int=0):
    """Parses and validates Job settings, then adds defaults and normalizes names."""
    job_data    = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    check_job_data(job_data, debug=max(debug-1,0))
    job_data    = fill_job_data(job_data)
    return job_data

def set_qc_data(content, section="&qc_data", isfile: bool=True, debug: int=0):
    """Parses and validates QC settings, then adds software-specific defaults."""
    qc_data     = Input_data(content=content, section=section, isfile=isfile, debug=debug)
    check_qc_data(qc_data, debug=max(debug-1,0))
    qc_data     = fill_qc_data(qc_data)
    return qc_data

def set_input_data(content, isfile: bool=True, debug: int=0):
    environment = set_environment_data(content, section="&environment", isfile=isfile, debug=debug)
    options     = set_options_data(content, section="&options", isfile=isfile, debug=debug)
    job_data    = set_job_data(content, section="&job_data", isfile=isfile, debug=debug)
    qc_data     = set_qc_data(content, section="&qc_data", isfile=isfile, debug=debug)
    return harmonize_input_data(environment, options, job_data, qc_data, debug=debug)

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

##############
def interpret_workflow_from_user(sys: object, job_data: object, debug: int=0):
    """Expands the user selection in job_data.workflow into a unique list of workflow names.

    workflow accepts a name/pattern or a list, matched case-insensitively against
    System sources. Shell-style *, ? and [...] patterns are supported; '*' selects
    all sources. Preserves selector order and source order within each match.
    Raises ValueError for invalid selectors or any selector matching no source.
    Returns the expanded list without changing job_data; branch membership is irrelevant.
    """
    if not hasattr(job_data,"workflow"):
        raise ValueError("INTERPRET_WORKFLOW_FROM_USER: job_data is missing 'workflow' input variable")
    if not hasattr(sys,"sources"):
        raise ValueError("INTERPRET_WORKFLOW_FROM_USER: system is missing 'sources'")

    selectors = [job_data.workflow] if isinstance(job_data.workflow,str) else job_data.workflow
    if not isinstance(selectors,list): 
        raise ValueError("INTERPRET_WORKFLOW_FROM_USER: job_data.workflow must be a source-name selector, or a list of selectors")

    names     = [sou.name.lower() for sou in sys.sources]
    workflows = []
    for selector in selectors:
        if not isinstance(selector,str) or not selector: 
            raise ValueError("INTERPRET_WORKFLOW_FROM_USER: job_data.workflow selectors must be nonempty strings")
        matches = [name for name in names if fnmatchcase(name,selector.lower())]
        if not matches: raise ValueError(f"INTERPRET_WORKFLOW_FROM_USER: No source matches workflow selector {selector!r}")
        for name in matches:
            if name not in workflows: workflows.append(name)

    if debug > 0: print(f"INTERPRET_WORKFLOW_FROM_USER: Selected {len(workflows)} workflows: {workflows}")
    return workflows
