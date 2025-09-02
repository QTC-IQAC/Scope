#################################
####  Contains the Cif Class ####
#################################
from .Parse_Cif import get_cif_diffraction_data, get_cif_authors, get_cif_journal 

###########
### CIF ###
###########
class cif(object):
    def __init__(self, name: str, path: str) -> None:
        self.type              = "cif" 
        self.version           = "1.0" 
        self.origin            = "created"
        self.name              = name
        self.path              = path
        self.read_bibliographic_data()

    ######
    def __repr__(self) -> None:
        to_print  = f'---------------------------------------------------\n'
        to_print += f'                   SCOPE .Cif file                 \n'
        to_print += f'---------------------------------------------------\n'
        to_print += f' Name                  = {self.name}\n'
        to_print += f' Path                  = {self.path}\n'
        if hasattr(self,"diff_temp"):        to_print += f' Diffraction Temp      = {self.diff_temp}\n'         
        if hasattr(self,"authors"):          to_print += f' Authors               = {self.authors}\n'           
        if hasattr(self,"journal_year"):     to_print += f' Year of Publication   = {self.journal_year}\n'      
        if hasattr(self,"journal_name"):     to_print += f' Journal Name          = {self.journal_name}\n'      
        if hasattr(self,"journal_volume"):   to_print += f' Journal Volume        = {self.journal_volume}\n'    
        if hasattr(self,"journal_page"):     to_print += f' Journal Page          = {self.journal_page}\n'      
        if hasattr(self,"cell"):             to_print += f' Has Associated Cell   = YES\n'
        return to_print

    ######
    def associate_cell(self, cell: object) -> None:
        self.cell      = cell
        return self.cell

    #######
    #def reset_path(self, new_path) -> None:
    #    if os.path.isdir(new_path): self.path = new_path

    ######
    def save(self, filepath: str=None):
        from Scope.Read_Write import save_binary
        if filepath is None: filepath = self.path
        save_binary(self, filepath)

    #######
    #def load_cell(self, cell_path: str) -> None:
    #    from Scope.Classes_Molecule import import_cell
    #    if os.path.isfile(cell_path): 
    #        try:
    #            self.cell      = import_cell(load_binary(self.cell_path))
    #            self.cell_path = cell_path
    #            return self.cell
    #        except:
    #            print(f"LOAD_CELL: Could not load cell from {cell_path}")  
    #            return None

    ##############################
    #### Bibliography options ####
    ##############################
    def read_bibliographic_data(self) -> None:
        self.diff_temp           = get_cif_diffraction_data(self.path) 
        self.authors             = get_cif_authors(self.path)
        year, name, volume, page = get_cif_journal(self.path)
        self.journal_year        = year 
        self.journal_name        = name 
        self.journal_volume      = volume 
        self.journal_page        = page 

    def get_search(self, verbose: bool=False, download: bool=True, debug: int=0):
        self.search = get_search(self, verbose=verbose, download=download, debug=debug) 
        return self.search

    def get_abstract(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.abstract = get_abstract(self.search, debug=debug)
        return self.abstract

    def get_title(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.title = get_title(self.search, debug=debug) 
        return self.title

    def get_doi(self, debug: int=0):
        if not hasattr(self,"search"): self.get_search(debug=debug)
        self.doi = get_doi(self.search, debug=debug) 
        return self.doi