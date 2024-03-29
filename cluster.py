## Cluster of immunoglobulin variable sequence using Hamming distance
import scipy
import scipy.cluster.hierarchy as sch
import matplotlib.pylab as plt
import operator
from itertools import imap
import numpy as np
import sys

##########################################################################
# Calculate the hamming distance between CDR3 pairs
def dist(str1, str2):
# Arguments:
#    strl: CDR3 amino sequence of sequence 1
#    str2: CDR3 amino sequence of sequence 2
    ne = operator.ne
    hamming = sum(imap(ne, str1, str2))
    return hamming

##########################################################################
# Construct the distance matric	
def pdist(X, metric):
# Arguments:
#    X: 
#    metric: methods for distance calculation
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2, ), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm

##########################################################################
# Hierarchical clustering
def single(sequence, key, cluster):
# Arguments:
#    sequence: dictionary of sequence
#    key；
#    cluster: dictionary of cluster
    y= pdist(sequence[key], dist)
    Z= sch.linkage(y, method='single')
    T = sch.fcluster(Z, 4, 'distance')
    for i in range(len(T)):
    if key+"-"+str(T[i]) in cluster:
        cluster[key+"-"+str(T[i])].append(sequence[key][i])
    else:
        cluster[key+"-"+str(T[i])]=[sequence[key][i]]

	

if __name__=='__main__':
    #___________________________________________________________________________
    # Combine the sequence with similar V gene, J gene and length of CDR3	
    with open (sys.argv[1], "r") as f:
        aa=f.readlines()
    
    cdrseq={}
    for i in range(len(aa)):
        temp=aa[i].split("\t")
        key_name=temp[1]+temp[4]+str(len(temp[0]))
        if (key_name) in cdrseq:
            cdrseq[key_name].append(aa[i])
        else:
	    cdrseq[key_name]=[aa[i]]
    
    # Cluster analysis
    seqcluster={}
    for key in cdrseq.keys():
        if len(cdrseq[key])==1:
	    seqcluster[key]=[cdrseq[key]]
        else:
            cluster={}
	    single(cdrseq, key, cluster)
	    seqcluster.update(cluster.items())
    
