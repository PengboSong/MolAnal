def new_cluster(atomn):
	cluster = []
	for i in range(atomn):
		cluster.append({"conf":i, "clust":i})
	return cluster

def link_structure(atomn, cutoff, rmsdSequence):
	cluster = new_cluster(atomn)
	sthChangeStatus = True
	while (sthChangeStatus is True):
		sthChangeStatus = False
		for item in rmsdSequence:
			if item.get("dist") < cutoff:
				xClust = cluster[item.get("x")].get("clust")
				yClust = cluster[item.get("y")].get("clust")
				diff = xClust - yClust
				if (diff != 0):
					sthChangeStatus = True
					if (diff > 0):
						cluster[item.get("y")].update({"clust":xClust})
					else:
						cluster[item.get("x")].update({"clust":yClust})
			else:
				break
	return cluster

def renum_cluster(cluster):
	clustID = 1
	for i in range(1, len(cluster)):
		if (cluster[i].get("clust") != cluster[i - 1].get("clust")):
			cluster[i - 1].update({"clust":clustID})
			clustID += 1
		else:
			cluster[i - 1].update({"clust":clustID})
	cluster[-1].update({"clust":clustID})

	clustList = [[] for i in range(clustID)]
	for item in cluster:
		clustList[item.get("clust") - 1].append(item.get("conf"))

	return cluster, clustList
