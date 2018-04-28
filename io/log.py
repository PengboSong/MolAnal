import os

def readlog(workdir = ".", startime = 0, timestep = 10):
	def check_selelog(s, length):
		try:
			int(s)
		except Exception:
			return False
		else:
			if 0 < int(s) <= length:
				return True
			else:
				return False
	# List all log file in the workdir
	logfiles = [f for f in os.listdir(workdir) if os.path.isfile(f) and os.path.splitext(f)[1]==".log"]
	for i in range(len(logfiles)):
		print('{0:>3}'.format(i+1)+": "+logfiles[i])
	# User input, select one id
	selelog = input("Enter selected log file id:")
	while not check_selelog(selelog, len(logfiles)):
		print("[Error] Invalid input option. Please enter again.")
		selelog = input("Enter selected log file id:")
	# Get log file name
	fname = os.path.splitext(logfiles[int(selelog)-1])[0]
	fpath = fname + ".log"
	csvpath = fname + ".csv"
	# Initialize data pack
	data = {}
	# Read log file
	with open(fpath, "r") as f:
		line = f.readline()
		while line[:6] != "frames":
			line = f.readline()
		ltitle = [x.strip() for x in line.strip().split("|")]
		for line in f.readlines():
			orgdist = [x.strip() for x in line.strip().split("|")]
			frame = int(orgdist[0])
			time = startime + timestep * frame
			dist = [float(x) for x in orgdist[1:]]
			data.update({time:dist})
	# Sort data by time order(increase)
	ltime = list(data.keys())
	ltime.sort()
	ldist = []
	for t in ltime:
		ldist.append(data.get(t))
	# Write csv
	with open(csvpath, "w") as f:
		f.write(";".join(ltitle) + '\n')
		f.writelines(str(ltime[i]) + ";" + ";".join(['{0:.3f}'.format(x) for x in ldist[i]]) + '\n' for i in range(len(ltime)))

def readxvg(workdir = "."):
	def check_selexvg(s, length):
		try:
			int(s)
		except Exception:
			return False
		else:
			if 0 < int(s) <= length:
				return True
			else:
				return False
	# List all xvg file in the workdir
	xvgfiles = [f for f in os.listdir(workdir) if os.path.isfile(f) and os.path.splitext(f)[1]==".xvg"]
	for i in range(len(xvgfiles)):
		print('{0:>3}'.format(i+1)+": "+xvgfiles[i])
	# User input, select one id
	selexvg = input("Enter selected xvg file id:")
	while not check_selexvg(selexvg, len(xvgfiles)):
		print("[Error] Invalid input option. Please enter again.")
		selexvg = input("Enter selected xvg file id:")
	# Get xvg file name
	fname = os.path.splitext(xvgfiles[int(selexvg)-1])[0]
	fpath = fname + ".xvg"
	csvpath = fname + ".csv"
	# Initialize data pack
	data = []
	# Read xvg file
	with open(fpath, "r") as f:
		line = f.readline()
		while line.startswith("@") or line.startswith("#"):
			line = f.readline()
		for line in f.readlines():
			data.append(line.strip().split())
	# Write csv
	with open(csvpath, "w") as f:
		f.writelines(";".join(line) + '\n' for line in data)
