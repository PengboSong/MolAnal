import os.path
import numpy as np
from custom.general import listFiles

def readlog(workdir = ".", start_time = 0.0, time_step = 10.0, write_new_table = True):
	data = []
	# Read log file
	fpath = listFiles(None, (".log", ".txt"), workdir)
	with open(fpath, "r") as f:
		ltitle = [x.strip() for x in f.readline().strip().split("|")]
		for line in f.readlines():
			orgdist = [x.strip() for x in line.strip().split("|")]
			frame_id = int(orgdist[0])
			real_time = start_time + time_step * frame_id
			line_data = [real_time]
			for v in orgdist[1:]:
				line_data.append(float(v))
			data.append(line_data)
	data = np.array(data)

	if write_new_table is True:
		csvpath = os.path.splitext(fpath)[0] + ".csv"
		# Write csv
		with open(csvpath, "w") as f:
			f.write(";".join(ltitle) + '\n')
			f.writelines(";".join(['{0:.3f}'.format(cell) for cell in row]) + '\n' for row in data)

	return data

def readxvg(workdir = ".", write_new_table = True):
	data = []
	# Read xvg file
	fpath = listFiles(None, ".xvg", workdir)
	with open(fpath, "r") as f:
		# Skip comment lines or format setting lines
		for line in f.readlines():
			if not line.startswith("@") and not line.startswith("#"):
				data.append([float(x) for x in line.strip().split()])
	data = np.array(data)

	if write_new_table is True:
		csvpath = os.path.splitext(fpath)[0] + ".csv"
		# Write csv
		with open(csvpath, "w") as f:
			f.writelines(";".join(['{0:.3f}'.format(cell) for cell in row]) + '\n' for row in data)
	return data
