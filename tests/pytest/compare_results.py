import json
import sys
import subprocess
import statistics

def readData(a, files):
    amin=[]
    amax=[]
    amin1=[]
    amax1=[]
    bad=[]
    ne1d=[]
    ne2d=[]
    ne3d=[]
    for f in files:
        for t in a[f]:
            if t['ne1d']>0:
                ne1d.append(t['ne1d'])
            if t['ne2d']>0:
                ne2d.append(t['ne2d'])
            if t['ne3d']>0:
                ne3d.append(t['ne3d'])
            if t['total_badness']>0.0:
                bad.append(t['total_badness'])
                if 'angles_tet' in t:
                    amin.append(t['angles_tet'][0])
                    amax.append(t['angles_tet'][1])
                if 'angles_trig' in t:
                    amin1.append(t['angles_trig'][0])
                    amax1.append(t['angles_trig'][1])
    return {
            "min tet angle":amin,
            "max tet angle" : amax,
            "min trig angle":amin1,
            "max trig angle" : amax1,
            "badness" : bad,
            "#edges" : ne1d,
            "#trigs" : ne2d,
            "#tets" : ne3d,
            }

import matplotlib.pyplot as plt

ref = 'master'
if len(sys.argv)>1:
    ref = sys.argv[1]

res = subprocess.run(['git','show','{}:./results.json'.format(ref)], capture_output=True)
s = json.loads(res.stdout.decode())

if len(sys.argv) > 2:
    ref2 = sys.argv[2]
    res = subprocess.run(['git','show','{}:./results.json'.format(ref2)], capture_output=True)
    s2 = res.stdout.decode()
else:
    ref2 = 'current'
    s2 = open('results.json','r').read()
s2 = json.loads(s2)

filenames = [f for f in s if f in s2]
data = readData(s, filenames)
data2 = readData(s2, filenames)

n = len(data)+1
fig,ax = plt.subplots(figsize=(10,7))
for i,d in enumerate(['min trig angle','min tet angle','max trig angle','max tet angle']):
    ax = plt.subplot(2,5,i+1)
    plt.title(d)
    ax.set_xticks([1,2])
    if len(data[d])==0 or len(data2[d])==0:
        continue

    plt.violinplot([data[d],data2[d]], showmedians=True)
    med = statistics.median(data[d])
    plt.hlines(med, 1,2, linestyle='dotted')
    if d=='badness':
        ax.set_yscale('log')
    ax.set_xticklabels([ref, ref2])


for i,d in enumerate(['badness','#edges','#trigs','#tets']):
    ax = plt.subplot(2,5,6+i)
    plt.title('difference '+d+' (in %)')
#     plt.violinplot([(y-x)/x for x,y in zip(data[d],data2[d])], showmedians=True)
    plt.boxplot([100.0*(y-x)/x for x,y in zip(data[d],data2[d])])
    plt.hlines(0.0, 0.5,1.5, linestyle='dotted')


# plt.savefig('comparison.png', dpi=100)
plt.show()
 
