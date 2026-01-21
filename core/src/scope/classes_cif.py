#################################
####  Contains the Cif Class ####
#################################
import sys
from scope.parse_general import search_string, read_lines_file 

###########
### CIF ###
###########
class Cif(object):
    def __init__(self, name: str, path: str) -> None:
        self.type              = "cif" 
        self.version           = "1.0" 
        self.origin            = "created"
        self.name              = name
        self.path              = path
        self.diff_temp         = get_cif_diffraction_data(path)
        self.sym_group         = get_symmetry_group(path)
        self.sym_ops           = get_symmetry_ops(path)
        self.get_biblio_data(path)

    ######
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print += f'                   SCOPE .Cif file                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Path                  = {self.path}\n'
        if self.sym_group is not None:       to_print += f' Symmetry Group        = {self.sym_group}\n'         
        if self.sym_ops is not None:         to_print += f' Symmetry Operations   = {self.sym_ops}\n'         
        if self.diff_temp is not None:       to_print += f' Diffraction Temp      = {self.diff_temp}\n'         
        if hasattr(self,"authors"):          to_print += f' Authors               = {self.authors}\n'           
        if hasattr(self,"journal_year"):     to_print += f' Year of Publication   = {self.journal_year}\n'      
        if hasattr(self,"journal_name"):     to_print += f' Journal Name          = {self.journal_name}\n'      
        if hasattr(self,"journal_volume"):   to_print += f' Journal Volume        = {self.journal_volume}\n'    
        if hasattr(self,"journal_page"):     to_print += f' Journal Page          = {self.journal_page}\n'      
        if hasattr(self,"cell"):             to_print += f' Has Associated Cell   = YES\n'
        else:                                to_print += f' Has Associated Cell   = NO\n'
        return to_print

    ######
    def associate_cell(self, cell: object) -> None:
        self.cell      = cell
        return self.cell

    ######
    def get_biblio_data(self, cifpath: str) -> None:
        self.authors   = get_cif_authors(cifpath)
        self.journal_year, self.journal_name, self.journal_volume, self.journal_page  = get_cif_journal(cifpath)

    ######
    def save(self, filepath: str=None):
        from scope.read_write import save_binary
        if filepath is None: filepath = self.path
        save_binary(self, filepath)

#############################
## Functions to Parse Cifs ##
#############################
def get_symmetry_group(cifpath: str):
    lines = read_lines_file(cifpath)
    sym_group_line, found       = search_string("_symmetry_space_group_name", lines, typ='first')
    if found: 
        assert len(lines[sym_group_line].split("'")) == 3, f"wrong block length: {lines[sym_group_line]}"
        sym_group = lines[sym_group_line].split("'")[1].rstrip()
    else:     
        print("Couldn't find symmetry group in cif:")
        sym_group = None
    return sym_group

def get_symmetry_ops(cifpath: str):
    lines = read_lines_file(cifpath)
    sym_ops_start, found1      = search_string("_symmetry_equiv_pos_as_xyz", lines, typ='first')
    sym_ops_end, found2        = search_string("_cell_length_a", lines, typ='first')
    assert sym_ops_end > sym_ops_start
    if found1 and found2: 
        sym_ops = []
        for l in lines[sym_ops_start+1:sym_ops_end]:
            assert len(l.split(" ")) == 2, f"wrong block length: {l}"
            sym_ops.append(l.split(" ")[1].rstrip())
    else:     
        print("Couldn't find symmetry operations in cif:")
        sym_ops = None
    return sym_ops

def get_cif_diffraction_data(cifpath: str):
    diff_temp = " "
    lines = read_lines_file(cifpath)
    diff_temp_line, found       = search_string("_diffrn_ambient_temperature", lines, typ='first')
    if found: 
        diff_temp = lines[diff_temp_line].split(" ")[1].rstrip()
    else:               
        print("Couldn't find diffraction temperature in cif:")
        diff_temp = None
    return diff_temp

def get_cif_authors(cifpath: str):
    lines = read_lines_file(cifpath)
    authors = " "
    authors_start, found1       = search_string("_publ_author_name", lines, typ='first')
    authors_end, found2         = search_string("_chemical_name_systematic", lines, typ='first')
    authors = []
    if found1 and found2:
        for i in range(authors_start+1, authors_end):
            aut = lines[i].rstrip().strip('"')
            authors.append(aut)
    else: print("Couldn't find authors in cif:", fil)
    return authors
 
def get_cif_journal(cifpath: str):
    lines = read_lines_file(cifpath)
    journal_year = journal_name = journal_volume = journal_page = " "
    journal_year_line, found3   = search_string("_journal_year", lines, typ='first')
    journal_name_line, found4   = search_string("_journal_name_full", lines, typ='first')
    journal_volume_line, found5 = search_string("_journal_volume", lines, typ='first')
    journal_page_line, found6   = search_string("_journal_page_first", lines, typ='first')                        
    if found3: 
        try:
            journal_year = lines[journal_year_line].split(" ")[1].rstrip()
        except Exception as exc: 
            print("Exception Reading Journal Year:", exc)
            print("Line is:", lines[journal_year_line])
    else: journal_year = '-'
    if found4:
        try:
            journal_name = lines[journal_name_line].split("'")[1].rstrip()
        except Exception as exc: 
            print("Exception Reading Journal Name:", exc)
            print("Line is:", lines[journal_name_line])
    else: journal_name = '-'
    if found5: 
        try:
            journal_volume = lines[journal_volume_line].split()[1].rstrip()
        except Exception as exc: 
            print("Exception Reading Journal Volume:", exc)
            print("Line is:", lines[journal_volume_line])
    else: journal_volume = '-'
    if found6: 
        try:
            journal_page = lines[journal_page_line].split()[1].rstrip()
        except Exception as exc: 
            print("Exception Reading Journal Page:", exc)
            print("Line is:", lines[journal_page_line])
    else: journal_page = '-'
    return journal_year, journal_name, journal_volume, journal_page

#########################
#### Other Functions ####
#########################
#def get_name_from_cif(cifpath: str):
#    lines = read_lines_file(cifpath)
#    journal_common, found   = search_string("_chemical_name_common",lines,type='first')
#    if int(journal_common) != 0: iscommon = True
#    else:                        iscommon = False
#    journal_chemname, found = search_string("_chemical_name_systematic",lines,type='first')
#    if iscommon:
#        chemname_start = int(journal_chemname+2)
#        chemname_end   = int(journal_common-2)
#    else:
#        journal_volume, found = search_string("_cell_volume",lines,type='first')
#        chemname_start = int(journal_chemname+2)
#        chemname_end   = int(journal_common-2)
#
#def get_volume_from_cif(cifpath: str):
#    lines = read_lines_file(cifpath)
#    journal_chemname, found = search_string("_cell_volume",lines,type='first')
