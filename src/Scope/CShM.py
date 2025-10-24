import os
import yaml
import numpy as np

# Load the ideal structures from YAML and check the structure
def get_cshm_ref(ref: str='OC-6', yaml_path: str=None, debug: int=0):
    if yaml_path is None:
        yaml_path = os.path.join(os.path.dirname(__file__), 'CShM.yaml')
    if not os.path.exists(yaml_path):
        print(f"Error: File '{yaml_path}' does not exist.")
        return None
    with open(yaml_path, 'r') as f:
        refs = yaml.safe_load(f)
    if ref not in refs: print(f"LOAD_REF: unknown {ref=}"); return None
    else:                                                   return np.array(refs[ref]) 
