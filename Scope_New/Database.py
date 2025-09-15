
def authors_to_references(aut_list):
    aut_str = ''
    if len(aut_list) == 1: aut_str += aut_list[0]
    elif len(aut_list) > 1:
        for idx in range(len(aut_list)-1):
            aut_str += aut_list[idx]+', '
        aut_str += 'and '
        aut_str += aut_list[-1]
    return aut_str
