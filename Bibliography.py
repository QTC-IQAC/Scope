import pybliometrics
from pybliometrics.scopus import ScopusSearch
from pybliometrics.scopus import AbstractRetrieval
import pandas as pd

pybliometrics.scopus.init()

def get_abstract(crys, debug: int=0):
    if debug > 0: print("GET_ABSTRACT: crys info:", crys.authors[0], crys.authors[-1], crys.journal_year, crys.journal_name, crys.journal_volume, crys.journal_page)
    query = f"AUTH ( {crys.authors[0].split('.')[-1]} ) AND AUTH ( {crys.authors[-1].split('.')[-1]} ) AND PAGEFIRST ( {crys.journal_page} ) AND PUBYEAR IS {crys.journal_year}"
    if debug > 0: print(query)

    s = ScopusSearch(query, verbose=True, download=True) 
    eids = s.get_eids()
    for e in eids:
        ab = AbstractRetrieval(e, view='FULL')
        print(ab.abstract)
        return ab.abstract

