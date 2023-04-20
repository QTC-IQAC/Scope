import os
from ast import literal_eval

#######################
##### INPUT CLASS #####
#######################
class QC_Input(object):
    def __init__(self, name: str, f_name = None):
        self.name = name
        if f_name is not None and isinstance(f_name, str):
            self.read(f_name)
 
    def read(self, f_name):
        dct = dict()
        with open(f_name, 'r') as f:
            for line in f.readlines():
                if line[0] == '#' or line == '\n': continue
                line = line.strip()
                data        = line.split('#')[0]
                key, value  = data.split(':')

                key         = key.strip()
                value       = value.strip()

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
        to_print += '---------------------------------------\n'
        to_print += string.format('Key', 'Data Type', 'Value')
        to_print += '---------------------------------------\n'
        for key in self.dct.keys():
            val = self.dct[key]
            to_print += string.format(key, str(type(val)), str(val))
        return to_print

