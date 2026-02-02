import numpy as np
from ast import literal_eval
from scope import constants

##################
### COLLECTION ###
##################
class Collection(object):
    def __init__(self, key: str, variable: str):
        self.type           = "collection"
        self.key            = key
        self.variable       = variable
        self.datas          = []

    def get_values(self):
        return np.asarray([data.value for data in self.datas])

    def get_variables(self):
        return np.asarray([getattr(d,self.variable) for d in self.datas])

    def add_data(self, data: object):
        if not hasattr(self,"units"):    self.units    = data.units
        if not hasattr(self,"function"): self.function = data.function
        if hasattr(self,"units") and self.units == data.units: 
            data._collection = self
            self.datas.append(data)

    def find_value_with_property(self, condition_name: str, condition_value):
        for idx, data in enumerate(self.datas):
            if hasattr(data, condition_name.lower()):
                value = getattr(data, condition_name.lower())
                if value == condition_value: return data

    def find_min(self):
        return self.datas[np.argmin(self.get_values())]

    def find_max(self):
        return self.datas[np.argmax(self.get_values())]

    def plot(self):
        import matplotlib.pyplot as plt
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(6, 4))
        x = self.get_variables()
        y = self.get_values()
        ax.scatter(x, y, color='blue', s=20, alpha=0.8, label='Data points')
        ax.set_xlabel(str(self.variable), fontsize=12)
        ax.set_ylabel(f"{self.key} (in {str(self.datas[0].units)})", fontsize=12)
        ax.set_title(f"{self.variable} vs {self.key} (in {self.datas[0].units})", fontsize=14, weight='bold')
        ax.legend(frameon=True, shadow=False)
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

####################
## Dunder Methods ##
####################
    def __len__(self):
        return int(len(self.datas))

    def __add__(self, other):
        assert isinstance(other, type(self))
        assert self.variable.lower() == other.variable.lower() # checks that they have the same variable 
        assert len(self) == len(other)                         # and length
        new_col = Collection(f"sum_{self.key}", self.variable)
        for idx, data1 in enumerate(self.datas):
            for jdx, data2 in enumerate(other.datas):
                prop1 = getattr(data1,self.variable.lower())
                prop2 = getattr(data2,self.variable.lower())
                if data1.units == data2.units and prop1 == prop2:
                    key   = f"sum_{self.key}"
                    value = data1.value + data2.value
                    units = data1.units
                    function = "collection.__sub__()"
                    new_data = Data(key, value, units, function)
                    new_data.add_property(self.variable, prop1, overwrite=True)
                    new_col.add_data(new_data)
        return new_col

    def __sub__(self, other):
        assert isinstance(other, type(self))
        assert self.variable.lower() == other.variable.lower() # checks that they have the same variable 
        assert len(self) == len(other)                         # and length
        new_col = Collection(f"delta_{self.key}", self.variable)
        for idx, data1 in enumerate(self.datas):
            for jdx, data2 in enumerate(other.datas):
                prop1 = getattr(data1,self.variable.lower())
                prop2 = getattr(data2,self.variable.lower())
                if data1.units == data2.units and prop1 == prop2:
                    key   = f"delta_{self.key}"
                    value = data1.value - data2.value
                    units = data1.units
                    function = "collection.__sub__()"
                    new_data = Data(key, value, units, function)
                    new_data.add_property(self.variable, prop1, overwrite=True)
                    new_col.add_data(new_data)
        return new_col

    def __repr__(self):
        if hasattr(self,"units"): 
            if self.units == 'kj': units = 'kJ/mol'
        to_print  = f'---------------------------------------------------\n'
        to_print += f'   COLLECTION OF DATA   = {self.key}                \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Variable              = {self.variable}\n'
        to_print += f' #Entries              = {len(self.datas)}\n'
        if len(self.datas) > 0:      to_print += f' First                 = {self.datas[0].value}\n'
        if len(self.datas) > 0:      to_print += f' Last                  = {self.datas[-1].value}\n'
        if hasattr(self,"units"):    to_print += f' Units                 = {self.units}\n'
        if hasattr(self,"function"): to_print += f' Function              = {self.function}\n'
        return to_print

############
### DATA ###
############
class Data(object):
    def __init__(self, key: str, value, units: str, function: str="Unknown", notes=None, debug: str=0):
        self.type          = "data"
        self.key           = key
        try: self.value    = literal_eval(value)
        except: self.value = value
        self.units         = units
        self.function      = function
        self.notes         = notes

    def add_property(self, name: str, value, overwrite: bool=False):
        if not hasattr(self, name):               setattr(self, name, value)
        elif   hasattr(self, name) and overwrite: setattr(self, name, value)
        else:  print("DATA.add_property: property already exists")

    def format(self):
        if self.units == 'kj': units = 'kJ/mol'
        elif self.units == 'au': units = 'au/mol'
        if type(self.value) == float: self.formatted = str(self.key+": "+str(f"{self.value:12.8f}")+" "+self.units)
        elif type(self.value) == int: self.formatted = str(self.key+": "+str(self.value)+" "+self.units)
        elif type(self.value) == str: self.formatted = str(self.key+": "+self.value+" "+self.units)
        elif self.value is None:      self.formatted = str(self.key+": None")

    def convert_to_units(self, new_units: str):
        if   self.units == 'au' and (new_units.lower() == 'kj' or new_units.lower() == 'kj/mol'):
            self.value = self.value * constants.har2kJmol
        elif self.units == 'kj' and new_units.lower() == 'au':
            self.value = self.value / constants.har2kJmol
        elif self.units == 'ry' and new_units.lower() == 'au':
            self.value = self.value * constants.ry2har
        elif self.units == 'au' and new_units.lower() == 'ry':
            self.value = self.value / constants.ry2har
        elif self.units == 'au' and new_units.lower() == 'ev':
            self.value = self.value * constants.har2eV
        elif self.units == 'ev' and new_units.lower() == 'au':
            self.value = self.value / constants.har2eV
        elif self.units == 'au' and new_units.lower() == 'cm':
            self.value = self.value * constants.har2cm
        elif self.units == 'cm' and new_units.lower() == 'au':
            self.value = self.value / constants.har2cm
        self.units = new_units
        self.format()
        return self

    def print_in_units(self, new_units: str):
        if new_units.lower() != self.units:
            if   self.units == 'au' and (new_units.lower() == 'kj' or new_units.lower() == 'kj/mol'):
                return f"{self.key}: {self.value * constants.har2kJmol:12.8f} {new_units}"
            elif self.units == 'kj' and new_units.lower() == 'au':
                return f"{self.key}: {self.value / constants.har2kJmol:12.8f} {new_units}"
            elif self.units == 'ry' and new_units.lower() == 'au':
                return f"{self.key}: {self.value * constants.r2har:12.8f} {new_units}"
            elif self.units == 'au' and new_units.lower() == 'ry':
                return f"{self.key}: {self.value / constants.r2har:12.8f} {new_units}"
            elif self.units == 'au' and new_units.lower() == 'ev':
                return f"{self.key}: {self.value * constants.har2eV:12.8f} {new_units}"
            elif self.units == 'ev' and new_units.lower() == 'au':
                return f"{self.key}: {self.value / constants.har2eV:12.8f} {new_units}"
            elif self.units == 'au' and new_units.lower() == 'cm':
                return f"{self.key}: {self.value * constants.har2cm:12.8f} {new_units}"
            elif self.units == 'cm' and new_units.lower() == 'au':
                return f"{self.key}: {self.value / constants.har2cm:12.8f} {new_units}"
        else:
            return str(self)
         
    def __repr__(self) -> None:
        if not hasattr(self,"formatted"): self.format()
        to_print =  f'{self.formatted}'
        return to_print

##################
### OPERATIONS ###
##################
def substract_collections(name: str, col1: object, col2: object, prop=None):
    ## Another substraction function, complementary of the dunder method
    assert col1.variable.lower() == col2.variable.lower()
    new_col = Collection(name, col1.variable)
    for idx, data1 in enumerate(col1.datas):
        for jdx, data2 in enumerate(col2.datas):
            if prop is not None: 
                if type(prop) != str: print("Substract_Collections: prop must be a string if not None"); return None
                prop1 = getattr(data1,prop)
                prop2 = getattr(data2,prop)
                if data1.units == data2.units and prop1 == prop2:
                    key = name
                    value = data1.value - data2.value
                    units = data1.units
                    function = "scope.Classes_Data.substract_collections()" 
                    new_data = Data(key, value, units, function)
                    new_data.add_property(prop, prop1, overwrite=True)
                    new_col.add_data(new_data)
    return new_col                
  

