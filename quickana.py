import os, time
import numpy as np

from gmx.structure.coord_matrix import load_coord_matrix
from gmx.structure.mol_matrix import load_mol_matrix, write_mol_matrix

from gmx.math.common import cluster, convertSeconds
from gmx.math.hbond import hbondmolgrp

from custom.progressbar import ProgressBar

def hbondgrp(moltypeA, moltypeB, start, end, lpdir = False, load_path = "frames", save_path = ".", prefix = "sysframe", suffix = "pdb"):
	# Initialize
	traj_prof = {}
	log = []

	frames_dir = os.path.join(save_path, "mols_%s-%s" %(moltypeA.lower(), moltypeB.lower()))
	if not os.path.isdir(frames_dir):
		os.mkdir(frames_dir)

	# Title
	title = '{0:>14}'.format("donor") + " | " + '{0:>14}'.format("hydro") + " | " + '{0:>14}'.format("acceptor") + " | " + '{0:>9}'.format("distance") + " | " + '{0:>6}'.format("angle")

	# Start timing
	startTime = time.time()

	# Load trajectory profile from the first frame
	first_frame = load_mol_matrix(os.path.join(load_path, "%s%d.%s" % (prefix, start, suffix)))
	# Remove coordinate and velocity terms
	# Add atom id term for each atom, begin counting
	# Atom ids would be count from 1
	atomn = 1
	moln = 0
	molids = list(first_frame.keys())
	molids.sort()
	for molid in molids:
		mol = first_frame.get(molid)
		mol_atomn = len(mol.get("atom"))
		mol_prof = {
		"name":mol.get("name"),
		"atom":mol.get("atom"),
		"atomid":[atomn + i for i in range(mol_atomn)]
		}
		atomn += mol_atomn
		moln += 1
		traj_prof.update({molid: mol_prof})

	prof_grp = cluster(traj_prof)
	assert moltypeA in prof_grp.keys() and moltypeB in prof_grp.keys(), "Wrong parameters. Function hbondgrp expects at least two parameters, both of them should be three-letter tags of corresponding molecules."

	# atomn equals to total counting of atoms minus initial 1
	# moln equals to total molecule number
	atomn -= 1
	# Only load coordinates data of each frame (including the first one) into a matrix (N*3)
	manage_load = ProgressBar(end - start + 1)
	manage_load.start("[Info] Start reading files and analysis.")
	framen = start
	frame_hbondns, frame_hbondmols = {}, {}
	for i in range(start, end + 1):
		frame_coord = load_coord_matrix(os.path.join(load_path, "%s%d.%s" % (prefix, i, suffix)))
		hbondn, orglog, molids = hbondmolgrp(traj_prof, frame_coord, moltypeA, moltypeB, lpdir)
		frame_hbondn = sum(hbondn)
		frame_hbondns.update({framen: frame_hbondn})
		if frame_hbondn > 0:
			log.append("<frame %d>" % framen)
			log.append(title)
			log.extend(orglog)
			frame_hbondmols.update({framen: molids})

		# Reconstruct the molecular matrix
		mol_matrix = {}
		moln = 0
		for j in molids:
			moln += 1
			mol_prof = traj_prof.get(j)
			mol = {"name":mol_prof.get("name"), "atom":mol_prof.get("atom")}
			mol.update({"coordinate":np.vstack([frame_coord[l-1] for l in mol_prof.get("atomid")])})
			mol_matrix.update({moln: mol})
		write_mol_matrix(mol_matrix, os.path.join(save_path, frames_dir,
			"%s+%s-frame%d.pdb" % (moltypeA.lower(), moltypeB.lower(), i)
		))

		framen += 1
		manage_load.forward(1)
	# End Timing
	endTime = time.time()
	manage_load.end("[Info] Hygrogen bond identification normally ends.")
	print("[Info] Total time usage = %s." % convertSeconds(endTime - startTime))

	if log:
		with open(os.path.join(save_path,
		"hbond_%s+%s-%s.txt" % (moltypeA.lower(), moltypeB.lower(), time.strftime("%Y%m%d_%H%M"))
		), "w") as f:
			f.writelines(line + '\n' for line in log)

	hbondnlog = ['{0:>6}'.format("frames") + " | " + '{0:>6}'.format("hbondn")]
	for i in range(self.trajstartn, self.trajendn + 1):
		hbondnlog.append('{0:>6}'.format(i) + " | " + '{0:>6}'.format(frame_hbondns.get(i)))
	with open(os.path.join(save_path,
		"hbondsum_%s+%s-%s.txt" % (moltypeA.lower(), moltypeB.lower(), time.strftime("%Y%m%d_%H%M"))
		), "w") as f:
		f.writelines(line + '\n' for line in hbondnlog)
