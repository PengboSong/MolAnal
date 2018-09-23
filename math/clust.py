def newCluster(atomn):
    cluster = []
    for i in range(atomn):
        cluster.append({"conf":i, "clust":i})
    return cluster

def linkStructure(atomn, cutoff, rmsdSequence):
    cluster = newCluster(atomn)
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

def renumCluster(cluster):
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

def oldSingleLinkage():
    # Cluster
	clusters = [[0 + self.trajstartn]]
	for i in range(1, self.trajn):
		newGroupFlag = True
		for group in clusters:
			# rmsdMatrix index equals to frame index minus start frame index, start from 0
			if all([self.rmsdMatrix[i][g - self.trajstartn] < cutoff for g in group]):
				group.append(i + self.trajstartn)
				newGroupFlag = False
				break
		if newGroupFlag is True:
			clusters.append([i + self.trajstartn])

	titleLine = '{0:>6}'.format("Groupn") + " | " + '{0:>24}'.format("Group members (Framen)")
	log = [titleLine]
	groupn = 0
	for group in clusters:
		groupn += 1
		log.append('{0:>6}'.format(groupn) + " | " + ', '.join([str(g) for g in group]))
