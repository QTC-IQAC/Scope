import sys

from Test_V3.Parse_General import search_string, read_lines_file 

def get_cif_diffraction_data(cifpath: str):
    diff_temp = " "
    lines = read_lines_file(cifpath)
    diff_temp_line, found = search_string("_diffrn_ambient_temperature", lines, typ='first')
    if found: 
        diff_temp = lines[diff_temp_line].split(" ")[1].rstrip()
    else: print("Couldn't find diffraction temperature in cif:", fil)
    return diff_temp

def get_cif_authors(cifpath: str):
    lines = read_lines_file(cifpath)
    authors = " "
    authors_start, found1 = search_string("_publ_author_name", lines, typ='first')
    authors_end, found2 = search_string("_chemical_name_systematic", lines, typ='first')
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
    journal_year_line, found3 = search_string("_journal_year", lines, typ='first')
    journal_name_line, found4 = search_string("_journal_name_full", lines, typ='first')
    journal_volume_line, found5 = search_string("_journal_volume", lines, typ='first')
    journal_page_line, found6 = search_string("_journal_page_first", lines, typ='first')                        
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
