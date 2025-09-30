import pybliometrics
from   pybliometrics.scopus import ScopusSearch
from   pybliometrics.scopus import AbstractRetrieval

#pybliometrics.scopus.init()

def get_query(cif, qtype: int=0, debug: int=0):
    """
    Generates a bibliographic query string based on the provided cif information and query type.

    Args:
        cif: A cif-class object containing bibliographic attributes such as authors, journal year, journal name, journal volume, and journal page.
        qtype (int, optional): Determines the format of the query.
            - 0: Uses first and last author, page, and year.
            - 1: Uses all authors' last names, page, and year.
            - 2: Uses all authors' first names, page, and year.
        debug (int, optional): If greater than 0, prints debug information.

    Returns:
        str: The formatted query string for bibliographic search.
    """
    if debug > 0: print("GET_ABSTRACT: cif info:", cif.authors[0], cif.authors[-1], cif.journal_year, cif.journal_name, cif.journal_volume, cif.journal_page)

    ## Type 0
    if qtype == 0:
        query = ""
        query += f"AUTH ( {cif.authors[0].split('.')[-1]} )" 
        query += f"AND AUTH ( {cif.authors[-1].split('.')[-1]} )" 
        query += f"AND PAGEFIRST ( {cif.journal_page} ) "
        query += f"AND PUBYEAR IS {cif.journal_year}"

    ## Type 1
    elif qtype == 1:
        query = ""
        for aut in cif.authors:
            query += f"AUTH ( {aut.replace('.',' ').split(' ')[-1]} )"
            query += f" AND "
        query += f"PAGEFIRST ( {cif.journal_page} ) "
        query += f"AND PUBYEAR IS {cif.journal_year}"

    ## Type 2
    elif qtype == 2:
        query = ""
        for aut in cif.authors:
            query += f"AUTH ( {aut.replace('.',' ').split(' ')[0]} )"
            query += f" AND "
        query += f"PAGEFIRST ( {cif.journal_page} ) "
        query += f"AND PUBYEAR IS {cif.journal_year}"

    ## Refining
    if "Junior" in query: query = query.replace("Junior","Jr.")
    if "junior" in query: query = query.replace("junior","Jr.")
    query = query.replace(u'\xa0', u' ')
    if debug > 0: print(f"GET_QUERY {query=}")
    return query

def get_search(cif, verbose: bool=False, download: bool=True, debug: int=0):
    for i in range(2):
        query = get_query(cif, qtype=i, debug=debug)
        search = ScopusSearch(query, verbose=verbose, download=download, debug=debug)
        size = search.get_results_size()
        if size == 1: return search 
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

def authors_to_references(aut_list):
    aut_str = ''
    if len(aut_list) == 1: aut_str += aut_list[0]
    elif len(aut_list) > 1:
        for idx in range(len(aut_list)-1):
            aut_str += aut_list[idx]+', '
        aut_str += 'and '
        aut_str += aut_list[-1]
    return aut_str
