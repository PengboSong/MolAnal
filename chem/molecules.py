# coding=utf-8

def hbondlist(moltype):
    """Get atom indexes for donor, hydrogen and acceptor identified by 3-char
    of molecules.

    For a classical hydrogen bond mode, a hydrogen can be recognized by
    donor-hydrogen/acceptor, showing that the heavy atom bonded with hydrogen
    is called donor and the other one is called acceptor.
    Attention: Atom indexes should start from 0, NOT 1.

    Returns:
            If 3-char code of molecule can be recognized, returns two lists.
            The first list containing tuples of (donor id, hydrogen id).
            The second one listing acceptor ids.
            If 3-char code can not be found, returns two empty lists.
    """
    hbond_terms = {
        "BCD": (
            [(0, 1), (4, 5), (10, 11), (19, 20), (22, 23), (25, 26), (33, 34), (36, 37), (39, 40), (47, 48), (50, 51),
             (53, 54), (61, 62), (64, 65), (67, 68), (75, 76), (78, 79), (81, 82), (87, 88), (90, 91), (96, 97)],
            [0, 4, 7, 10, 12, 14, 17, 19, 22, 25, 28, 31, 33, 36, 39, 42, 45, 47,
                50, 53, 56, 59, 61, 64, 67, 70, 73, 75, 78, 81, 84, 87, 90, 93, 96]
        ),
        "C8A": ([(2, 11)], [0, 2]),
        "C8O": ([(8, 9)], [8]),
        "C8N": ([(8, 9), (8, 10)], [8]),
        "NO3": ([], [1, 2, 3]),
        "DCF": ([(8, 22)], [18, 19]),
        "LAS": ([], [7, 8, 9]),
        "SOL": ([(0, 1), (0, 2)], [0]),
    }
    if moltype in hbond_terms:
        return hbond_terms.get(moltype)
    else:
        return ([], [])


def hbondlist_direction(moltype):
    """Get detailed information to further identify a hydrogen bond, including
    donor atom hybridization, donor atom index and its surrounding atoms.

    Assuming that hydrogen bond can only be formed around axes determined by
    acceptor hybridization. For each donor atoms, hybridization, atom index and
    its surrounding atoms will be packed as follows:
    HYBRIDIZATION DONOR_ATOM_ID CONNECT_ATOM1_ID CONNECT_ATOM2_ID ...
    For example, for a donor atom 0 connected by atom 1, which should be a "sp2"
    hybridization, it should be packed like this: ("sp2", 0, 1)
    Attention: Atom indexes should start from 0, NOT 1.

    Returns:
            If 3-char code of molecule can be recognized, returns a list containing
            tuples of (hybridization, donor atom id, connect atom ids).
            If 3-char code can not be found, returns an empty list.
    """
    hbond_direction_terms = {
        "NO3": [("sp2", 1, 0), ("sp2", 2, 0), ("sp2", 3, 0)],
        "SOL": [("sp3", 0, 1, 2)],
    }
    if moltype in hbond_direction_terms:
        return hbond_direction_terms.get(moltype)
    else:
        return []


def conect(moltype):
    """Get PDB-like connect terms identified by 3-char of molecules.

    PDB format file records bond/connect information by terms as follows:
    CENTRAL_ATOM_ID CONNECT_ATOM1_ID CONNECT_ATOM2_ID CONNECT_ATOM3_ID ...
    Connect information is packed into a tuple in the same order.	
    Attention: Atom indexes should start from 1.

    Returns:
            If 3-char code of molecule can be recognized, returns a list containing
            tuples with connect atom indexes packed.
            If 3-char code can not be found, returns an empty list.
    """
    conect_terms = {
        "QAC": [(1, 2, 6, 10, 14), (2, 3, 4, 5, 1), (3, 2), (4, 2), (5, 2), (6, 7, 8, 9, 1), (7, 6), (8, 6), (9, 6), (10, 11, 12, 13, 1), (11, 10), (12, 10), (13, 10), (14, 15, 18, 19, 1), (15, 14, 16, 20, 21), (16, 15, 17, 22, 23), (17, 16, 24, 25, 26), (18, 14), (19, 14), (20, 15), (21, 15), (22, 16), (23, 16), (24, 17), (25, 17), (26, 27, 30, 31, 17), (27, 26, 28, 32, 33), (28, 27, 29, 34, 35), (29, 28, 36, 37, 38), (30, 26), (31, 26), (32, 27), (33, 27), (34, 28), (35, 28), (36, 29), (37, 29), (38, 39, 42, 43, 29), (39, 38, 40, 44, 45), (40, 39, 41, 46, 47), (41, 40, 48, 49, 50), (42, 38), (43, 38), (44, 39), (45, 39), (46, 40), (47, 40), (48, 41), (49, 41), (50, 51, 54, 55, 41), (51, 50, 52, 56, 57), (52, 51, 53, 58, 59), (53, 52, 60, 61, 63), (54, 50), (55, 50), (56, 51), (57, 51), (58, 52), (59, 52), (60, 53), (61, 53), (62, 63, 64, 65, 66), (63, 62, 67, 68, 53), (64, 62), (65, 62), (66, 62), (67, 63), (68, 63)],
        "NO3": [(1, 2, 3, 4), (2, 1), (3, 2), (4, 2)],
    }
    if moltype in conect_terms:
        return conect_terms.get(moltype)
    else:
        return []
