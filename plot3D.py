from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

v=[];d=[];j=[];p=[]

total_data<-open("out.vdj.table.patient.fasta.txt","rb")
with open(total_data) as f:
    data=f.readlines()

for line in range(len(data)):
    v.append(line.split("\t")[0])
    d.append(line.split("\t")[1])
    j.append(line.split("\t")[2])
    f.append(line.split("\t")[3])

### Calculation of the frequency of each immunoglobulin gene family!
def count_list(l):
# Args:
# l
	d={}
	for k in l:
		if d.has_key(k):
			d[k]+=1
		else:
			d[k]=1
	return d

### transform the immunoglobulin gene family into numeric data!
def trans(family):
#    family=[]
    unique=[i for i in set(family)]
    sort_name=sorted(unique,reverse=False)
    h={}
    for i in range(0,len(sort_name)):
        h[sort_name[i]]=i		
    for i in range(0,len(family)):
	if family[i] in h:
		family[i]=h[family[i]]

### 3D plot
fig=plt.figure()
ax=fig.add_subplot(111,projection="3d")
ax.set_xlabel('V family')
ax.set_ylabel('D family')
ax.set_zlabel('J family')
plt.titlt("85-426 Heavy Chain")
plt.show()    


