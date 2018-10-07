def hbondlist(moltype):
	# Attention: atom ids in each list should start from 0, not 1
	# Add new molecular type here when trying to analysis a new system
	# Each value of the given molecular type has two terms: the former term is the donor list, with the donor atom be the first and the hydrogen be the last, and the latter one is the acceptor list
	hbondTerms = {
	"BCD": (
		[(0, 1), (4, 5), (10, 11), (19, 20), (22, 23), (25, 26), (33, 34), (36, 37), (39, 40), (47, 48), (50, 51), (53, 54), (61, 62), (64, 65), (67, 68), (75, 76), (78, 79), (81, 82), (87, 88), (90, 91), (96, 97)],
		[0, 4, 7, 10, 12, 14, 17, 19, 22, 25, 28, 31, 33, 36, 39, 42, 45, 47, 50, 53, 56, 59, 61, 64, 67, 70, 73, 75, 78, 81, 84, 87, 90, 93, 96]
	),
	"C8A": ([(2, 11)], [0, 2]),
	"C8O": ([(8, 9)], [8]),
	"C8N": ([(8, 9), (8, 10)], [8]),
	"QAC": ([], []),
	"NO3": ([], [1, 2, 3]),
	"DCF": ([(8, 22)], [18, 19]),
	"LAS": ([], [7, 8, 9]),
	"SOL": ([(0, 1), (0, 2)], [0])
	}
	if not (moltype in hbondTerms):
		raise ValueError("Molecular type %s has no corresponding data term.." % moltype)
	else:
		return hbondTerms.get(moltype)

def hbondlist_direction(moltype):
	# Attention: atom ids in each list should start from 0, not 1
	# Add new molecular type here when trying to analysis a new system
	hbondDirectionTerms = {
	"BCD": [],
	"C8A": [],
	"C8O": [],
	"C8N": [],
	"QAC": [],
	"NO3": [("sp2", 1, 0), ("sp2", 2, 0), ("sp2", 3, 0)],
	"DCF": [],
	"LAS": [],
	"SOL": [("sp3", 0, 1, 2)]
	}
	if not (moltype in hbondDirectionTerms):
		raise ValueError("Molecular type %s has no corresponding data term.." % moltype)
	else:
		return hbondDirectionTerms.get(moltype)
