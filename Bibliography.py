import pybliometrics
from pybliometrics.scopus import ScopusSearch
from pybliometrics.scopus import AbstractRetrieval
import pandas as pd

pybliometrics.scopus.init()

def get_query(crys, debug: int=0):
    if debug > 0: print("GET_ABSTRACT: crys info:", crys.authors[0], crys.authors[-1], crys.journal_year, crys.journal_name, crys.journal_volume, crys.journal_page)
    query = f"AUTH ( {crys.authors[0].split('.')[-1]} ) AND AUTH ( {crys.authors[-1].split('.')[-1]} ) AND PAGEFIRST ( {crys.journal_page} ) AND PUBYEAR IS {crys.journal_year}"
    if debug > 0: print(query)
    return query

def get_abstract(query: str, debug: int=0):
    s = ScopusSearch(query, verbose=True, download=False)
    eids = s.get_eids()
    for e in eids:
        ab = AbstractRetrieval(e, view='FULL')
        return ab.abstract

def get_doi(query: str, debug: int=0):
    s = ScopusSearch(query, verbose=True, download=False)
    eids = s.get_eids()
    for e in eids:
        ab = AbstractRetrieval(e, view='FULL')
        return ab.doi

def get_title(query: str, debug: int=0):
    s = ScopusSearch(query, verbose=True, download=False)
    eids = s.get_eids()
    for e in eids:
        ab = AbstractRetrieval(e, view='FULL')
        return ab.title
