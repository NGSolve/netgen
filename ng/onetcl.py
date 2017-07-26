import glob
import functools

# load all tcl files in current folder
tclfiles = {}
for fname in glob.glob('*.tcl'):
    tclfiles[fname] = open(fname,'r').read()

# do a topological sorting (such that if a.tcl is including b.tcl,
# a will come after b in the sorted list)
fnames = list(tclfiles.keys())
fnames.sort(key=functools.cmp_to_key( lambda x,y: tclfiles[x].find('/'+y) ))

# replace all occurrences of 'source bla.tcl' with the code of bla.tcl
for f in fnames:
    for g in fnames:
        if(tclfiles[f].find('/'+g) >= 0):
            tclfiles[f] = tclfiles[f].replace("source ${ngdir}/"+g, tclfiles[g])

# write a cpp file containing the result of ng.tcl
onetclcpp = open("onetcl.cpp",'w')
onetclcpp.write('const char * ngscript[] = {""'+'\n');

# make sure to remove comments (and if lines with comments end with '\' also the next line(s) )
skip_next = False # flag to indicate that the next line should be removed
for line in tclfiles["ng.tcl"].split('\n'):
    line = line.strip()
    if len(line)==0:
        skip_next = False
        continue
    if skip_next:
        # skip as long as lines end with '\'
        skip_next = line[-1]=='\\'
        continue
    if(line.find('#')>-1):
        # comment found (not necessarily the whole line)
        skip_next = line[-1]=='\\'
        line = line[:line.find('#')]
    if len(line)>0:
        s = ',"' + line.replace('\\', r'\\').replace('"', r'\"') + '\\n"\n'
        onetclcpp.write(s)

onetclcpp.write(', nullptr\n');
onetclcpp.write('};\n');
onetclcpp.close();
