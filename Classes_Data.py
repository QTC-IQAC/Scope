from ast import literal_eval

##################
### COLLECTION ###
##################
class collection(object):
    def __init__(self, key: str):
        self.type  = "collection"
        self.key   = key
        self.datas = []

    def get_values(self):
        values = []
        for data in self.datas:
            values.append(data.value)
        return values 
        #return list([data.value for data in self.datas])

    def add_data(self, data: object):
        if not hasattr(self,"units"):    self.units    = data.units
        if not hasattr(self,"function"): self.function = data.function
        if hasattr(self,"units") and self.units == data.units: 
            data._colection = self
            self.datas.append(data)

    def find_value_with_property(self, condition_name: str, condition_value):
        for idx, data in enumerate(self.datas):
            if hasattr(data, condition_name):
                value = getattr(data, condition_name)
                if value == condition_value: return data

    ### Should be changed to __repr__
    def format(self):
        if self.units == 'kj': units = 'kJ/mol'
        to_print  = f'---------------------------------------------------\n'
        to_print += f'   COLECTION OF DATA   = {self.key}                \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' _object.type          = {self._object.type}\n'
        to_print += f' _object.keyword       = {self._object.keyword}\n'
        to_print += f' #Entries              = {len(self.datas)}\n'
        if len(self.datas) > 0: to_print += f' First                 = {self.datas[0].value}\n'
        if len(self.datas) > 0: to_print += f' Last                  = {self.datas[-1].value}\n'
        if hasattr(self,"units"):    to_print += f' Units                 = {self.units}\n'
        if hasattr(self,"function"): to_print += f' Function              = {self.function}\n'
        print(to_print)

############
### DATA ###
############
class data(object):
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
        else:  print("add_property: property already exists")

    def format(self):
        if self.units == 'kj': units = 'kJ/mol'
        elif self.units == 'au': units = 'au/mol'
        if type(self.value) == float: self.formatted = str(self.key+": "+str(f"{self.value:12.8f}")+" "+self.units)
        elif type(self.value) == int: self.formatted = str(self.key+": "+str(self.value)+" "+self.units)
        elif type(self.value) == str: self.formatted = str(self.key+": "+self.value+" "+self.units)
        elif self.value is None:      self.formatted = str(self.key+": None")
         
    def __repr__(self) -> None:
        if not hasattr(self,"formatted"): self.format()
        to_print =  f'{self.formatted}\n'
        #to_print += f'------------------------------------------\n'
        #to_print += f'Extracted from function: {self.function}  \n'
        #to_print += f'Notes: {self.notes}                       \n'
        #to_print += f'------------------------------------------\n'
        return to_print

##################
### OPERATIONS ###
##################
def substract_collections(name: str, col1: object, col2: object, prop=None):
    new_col = collection(name)
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
                    function = "substract_collections" 
                    new_data = data(key, value, units, function)
                    new_data.add_property(prop, prop1, overwrite=True)
                    new_col.add_data(new_data)
    return new_col                
  

