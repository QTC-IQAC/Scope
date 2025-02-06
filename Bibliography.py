import pybliometrics
from pybliometrics.scopus import ScopusSearch
from pybliometrics.scopus import AbstractRetrieval
import pandas as pd

pybliometrics.scopus.init()

def get_query(crys, qtype: int=0, debug: int=0):
    if debug > 0: print("GET_ABSTRACT: crys info:", crys.authors[0], crys.authors[-1], crys.journal_year, crys.journal_name, crys.journal_volume, crys.journal_page)
    ## Type 0
    if qtype == 0:
        query = ""
        query += f"AUTH ( {crys.authors[0].split('.')[-1]} )" 
        query += f"AND AUTH ( {crys.authors[-1].split('.')[-1]} )" 
        query += f"AND PAGEFIRST ( {crys.journal_page} ) "
        query += f"AND PUBYEAR IS {crys.journal_year}"
    ## Type 1
    elif qtype == 1:
        query = ""
        for aut in crys.authors:
            query += f"AUTH ( {aut.replace('.',' ').split(' ')[-1]} )"
            query += f" AND "
        query += f"PAGEFIRST ( {crys.journal_page} ) "
        query += f"AND PUBYEAR IS {crys.journal_year}"
    ## Type 2
    elif qtype == 2:
        query = ""
        for aut in crys.authors:
            query += f"AUTH ( {aut.replace('.',' ').split(' ')[0]} )"
            query += f" AND "
        query += f"PAGEFIRST ( {crys.journal_page} ) "
        query += f"AND PUBYEAR IS {crys.journal_year}"

    ## Polish
    if "Junior" in query: query = query.replace("Junior","Jr.")
    if "junior" in query: query = query.replace("junior","Jr.")
    query = query.replace(u'\xa0', u' ')
    if debug > 0: print(f"GET_QUERY {query=}")
    return query

def get_search(crys, verbose: bool=False, download: bool=True, debug: int=0):
    for i in range(2):
        query = get_query(crys, qtype=i, debug=debug)
        search = ScopusSearch(query, verbose=verbose, download=download, debug=debug)
        size = search.get_results_size()
        if size == 1: return search 
        #if debug > 0: print(f"GET_SEARCH {search=} with {query=}")
    return None

def get_abstract(search, debug: int=0):
    if search is None: return None
    if search.get_results_size() == 0: return None
    for e in search.get_eids():
        ab = AbstractRetrieval(e, view='FULL')
        return ab.abstract

def get_doi(search, debug: int=0):
    if search is None: return None
    if search.get_results_size() == 0: return None
    for e in search.get_eids():
        ab = AbstractRetrieval(e, view='FULL')
        return ab.doi
    
def get_title(search, debug: int=0):
    if search is None: return None
    if search.get_results_size() == 0: return None
    for e in search.get_eids():
        ab = AbstractRetrieval(e, view='FULL')
        return ab.title
    
#def get_abstract(query: str, verbose: bool=False, download: bool=True, debug: int=0):
#    s = ScopusSearch(query, verbose=verbose, download=download)
#    eids = s.get_eids()
#    for e in eids:
#        ab = AbstractRetrieval(e, view='FULL')
#        return ab.abstract
#
#def get_doi(query: str, verbose: bool=False, download: bool=True, debug: int=0):
#    s = ScopusSearch(query, verbose=verbose, download=download)
#    eids = s.get_eids()
#    for e in eids:
#        ab = AbstractRetrieval(e, view='FULL')
#        return ab.doi
#
#def get_title(query: str, verbose: bool=False, download: bool=True, debug: int=0):
#    s = ScopusSearch(query, verbose=verbose, download=download)
#    eids = s.get_eids()
#    for e in eids:
#        ab = AbstractRetrieval(e, view='FULL')
#        return ab.title
