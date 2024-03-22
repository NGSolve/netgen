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
    file=[]
    for f in files:
        if f == 'cylinder.geo':
            continue
        for t in a[f]:
            if t['ne1d']>0:
                ne1d.append(t['ne1d'])
            if t['ne2d']>0:
                ne2d.append(t['ne2d'])
            if t['ne3d']>0:
                ne3d.append(t['ne3d'])
            if t['total_badness']>0.0:
                bad.append(t['total_badness'])
                file.append(f)
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
            "file" : file,
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

assert(len(data) == len(data2))

w = 90
GREEN = '\033[92m'
RED = '\033[91m'
RESET = '\033[0m'

for bad1,bad2, f1, f2 in zip(data['badness'], data2['badness'], data['file'], data2['file']):
    assert f1==f2

    diff = f"{100*(bad2-bad1)/bad1:+.2f}%"
    if bad2>0 and bad2>1.2*bad1:
        print(f"{RED}badness {f1} got worse:  {bad1} -> {bad2}".ljust(w) + diff + RESET)
    if bad2>0 and bad2<0.8*bad1:
        print(f"{GREEN}badness {f1} got better: {bad1} -> {bad2}".ljust(w) + diff + RESET)

for bad1,bad2, f1, f2 in zip(data['#trigs'], data2['#trigs'], data['file'], data2['file']):
    assert f1==f2
    diff = f"{100*(bad2-bad1)/bad1:+.2f}%"
    if bad2>0 and bad2>1.2*bad1:
        print(f"{RED}ntrigs {f1} got worse:  {bad1} -> {bad2}".ljust(w) + diff + RESET)
    if bad2>0 and bad2<0.8*bad1:
        print(f"{GREEN}ntrigs {f1} got better: {bad1} -> {bad2}".ljust(w) + diff + RESET)

n = len(data) + 1
fig, ax = plt.subplots(figsize=(15, 7))  # Adjust figsize as needed
plt.xticks([])
plt.yticks([])
ax.yaxis.grid(False)
ax.xaxis.grid(False)
for i, d in enumerate(['min trig angle', 'min tet angle', 'max trig angle', 'max tet angle']):
    plt.subplot(2, 4, i + 1)  # Remove ax = 
    plt.title(d)
    # plt.xticks([1, 2])
    if len(data[d]) == 0 or len(data2[d]) == 0:
        continue

    plt.violinplot([data[d], data2[d]], showmedians=True)
    med = statistics.median(data[d])
    plt.hlines(med, 1, 2, linestyle='dotted')
    if d == 'badness':
        plt.yscale('log')
    plt.xticks([1, 2], [ref, ref2])

for i, d in enumerate(['badness', '#edges', '#trigs', '#tets']):
    plt.xticks([])
    plt.subplot(2, 4, 5 + i)
    plt.title('difference ' + d + ' (in %)')
    plt.boxplot([100.0 * (y - x) / x for x, y in zip(data[d], data2[d])])
    plt.hlines(0.0, 0.5, 1.5, linestyle='dotted')

plt.tight_layout()  # Adjust layout


# plt.savefig('comparison.png', dpi=100)
plt.show()
 
