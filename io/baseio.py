# coding=utf-8


from enum import Enum
import re


class GMXDomains(Enum):
    # Domain identifier : Short name
    RECORD = "record"
    ATOMID = "atom_id"
    ATOMNM = "atom_name"
    ELEMENT = "element"
    RESID = "res_id"
    RESNM = "res_name"
    MOLID = "mol_id"
    MOLNM = "mol_name"
    X = "x"
    Y = "y"
    Z = "z"
    XYZ = "coordinate"
    VX = "vx"
    VY = "vy"
    VZ = "vz"
    VXYZ = "velocity"
    LOC = "location"
    CHAIN = "chain"
    CODE = "code"
    OCCUP = "occup"
    TEMP = "temp_factor"
    SEGMENT = "segment"
    OATOMID = "obsolete_atom_id"


DEF_VALUES = {
    "character": '',
    "integer": 0,
    "real": 0.,
}


def gmx_readline(line, rown, domains, selected_domains=None):
    '''Read one line that data are separated by columns.

    Extract useful values by column numbers.

    Args:
        line(str): One line from file.
        rown(int): Index of this line.
        domains(dict[GMXDomains : tuple(str, int, int, type, str, str)]):
            Details of domains recording long name, start/end location,
            right Python type, type name and alignment.
        selected_domains(list[GMXDomains]): List of selected domain indexes.

    Returns:
        A dict of extracted domains with GMXDomains as keys.

    Raises:
        ValueError: Domain contents can not be recognized.
        AssertionError: Received unparsed line or invalid row number.
    '''

    assert isinstance(line, str), "Can not parse line {}".format(line)
    assert isinstance(rown, int), "Invalid row number {} given".format(rown)

    if not selected_domains:
        selected_domains = domains.keys()

    domain_values = {}
    for domain in selected_domains:
        if domain in GMXDomains:
            lname, sloc, eloc, vtype, typename, _ = domains[domain]
            vcontent = line[sloc:eloc].strip()
            if len(vcontent) == 0:
                vcontent = DEF_VALUES[typename]
            try:
                value = vtype(vcontent)
                domain_values.update({domain: value})
            except Exception:
                raise ValueError("{0:s} at row {1:d} column {2:d}-{3:d} can not "
                                    "be recognized as {4:s}.".format(lname, rown, sloc, eloc, typename))

    return domain_values


def readline_sep(line, split_symbol=';', comment_symbol='#'):
    content = re.split(comment_symbol, line, maxsplit=1)[0].strip()
    if content:
        return re.split(split_symbol, content)
    else:
        return []


def gmx_readline_sep(line, rown, domains, selected_domains=None, split_symbol=';', comment_symbol='#'):
    '''Read one line that data are separated by separators.

    Extract useful values by blocks separated by separators.

    Args:
        line(str): One line from file.
        rown(int): Index of this line.
        domains(dict[GMXDomains : tuple(str, int, type, str, str)]):
            Details of domains recording long name, block location,
            right Python type, type name and alignment.
        indexes(list[GMXDomains]): List of selected domain indexes.
        split_symbol(str): Collection of separator symbols.
        comment_symbol(str): Collection of comment symbols after which are
                             line comments.

    Returns:
        A dict of extracted domains with GMXDomains as keys.

    Raises:
        ValueError: Domain contents can not be recognized.
        AssertionError: Received unparsed line or invalid row number.
    '''

    assert isinstance(line, str), "Can not parse line {}".format(line)
    assert isinstance(rown, int), "Invalid row number {} given".format(rown)

    if not selected_domains:
        selected_domains = domains.keys()

    lineblocks = readline_sep(
        line, split_symbol=split_symbol, comment_symbol=comment_symbol)

    domain_values = {}
    for domain in selected_domains:
        if domain in GMXDomains:
            lname, loc, vtype, typename, _ = domains[domain]
            try:
                value = vtype(lineblocks[loc].strip())
                domain_values.update({domain: value})
            except Exception:
                raise ValueError("{0:s} at row {1:d} column {2:d} can not "
                                 "be recognized as {3:s}.".format(lname, rown, loc, typename))

    return domain_values
