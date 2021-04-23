def hbondlist(moltype):
	# Attention: atom ids in each list should start from 0, not 1
	# Add new molecular type here when trying to analysis a new system
	# Each value of the given molecular type has two terms: the former term is the donor list, with the donor atom be the first and the hydrogen be the last, and the latter one is the acceptor list
	hbond_terms = {
	"BCD": (
		[(0, 1), (4, 5), (10, 11), (19, 20), (22, 23), (25, 26), (33, 34), (36, 37), (39, 40), (47, 48), (50, 51), (53, 54), (61, 62), (64, 65), (67, 68), (75, 76), (78, 79), (81, 82), (87, 88), (90, 91), (96, 97)],
		[0, 4, 7, 10, 12, 14, 17, 19, 22, 25, 28, 31, 33, 36, 39, 42, 45, 47, 50, 53, 56, 59, 61, 64, 67, 70, 73, 75, 78, 81, 84, 87, 90, 93, 96]
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
	# Attention: atom ids in each list should start from 0, not 1
	# Add new molecular type here when trying to analysis a new system
	hbond_direction_terms = {
	"NO3": [("sp2", 1, 0), ("sp2", 2, 0), ("sp2", 3, 0)],
	"SOL": [("sp3", 0, 1, 2)],
	}
	if moltype in hbond_direction_terms:
		return hbond_direction_terms.get(moltype)
	else:
		return []

def conect(moltype):
	# Attention: atom ids in each list should start from 1
	conect_terms = {
	"QAC": [(1, 2, 6, 10, 14), (2, 3, 4, 5, 1), (3, 2), (4, 2), (5, 2), (6, 7, 8, 9, 1), (7, 6), (8, 6), (9, 6), (10, 11, 12, 13, 1), (11, 10), (12, 10), (13, 10), (14, 15, 18, 19, 1), (15, 14, 16, 20, 21), (16, 15, 17, 22, 23), (17, 16, 24, 25, 26), (18, 14), (19, 14), (20, 15), (21, 15), (22, 16), (23, 16), (24, 17), (25, 17), (26, 27, 30, 31, 17), (27, 26, 28, 32, 33), (28, 27, 29, 34, 35), (29, 28, 36, 37, 38), (30, 26), (31, 26), (32, 27), (33, 27), (34, 28), (35, 28), (36, 29), (37, 29), (38, 39, 42, 43, 29), (39, 38, 40, 44, 45), (40, 39, 41, 46, 47), (41, 40, 48, 49, 50), (42, 38), (43, 38), (44, 39), (45, 39), (46, 40), (47, 40), (48, 41), (49, 41), (50, 51, 54, 55, 41), (51, 50, 52, 56, 57), (52, 51, 53, 58, 59), (53, 52, 60, 61, 63), (54, 50), (55, 50), (56, 51), (57, 51), (58, 52), (59, 52), (60, 53), (61, 53), (62, 63, 64, 65, 66), (63, 62, 67, 68, 53), (64, 62), (65, 62), (66, 62), (67, 63), (68, 63)],
	"NO3": [(1, 2, 3, 4), (2, 1), (3, 2), (4, 2)],
	}
	if moltype in conect_terms:
		return conect_terms.get(moltype)
	else:
		return []
