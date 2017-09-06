from netgen.meshing import *

def ReadGmsh (filename):
    f = open(filename, 'r')

    mesh = Mesh(dim=3)
    
    pointmap = { }
    facedescriptormap = { } 
                         
    
    while True:
        line = f.readline()
        if line == "": break
            
        if line.split()[0]=="$Nodes":
            num = int(f.readline().split()[0])
            print ("reading",num,"nodes")
            for i in range(num):
                line = f.readline()
                nodenum,x,y,z = line.split()[0:4]
                pnum = mesh.Add (MeshPoint(Pnt(float(x),float(y),float(z))))
                pointmap[int(nodenum)] = pnum

        if line.split()[0]=="$Elements":
            num = int(f.readline().split()[0])
            print ("reading",num,"elements")
            
            for i in range(num):
                line = f.readline().split()
                elmnum = int(line[0])
                elmtype = int(line[1])
                numtags = int(line[2])
                tag1 = int(line[3])
                    
                if elmtype == 2:   # 3-node trig
                    num_nodes = 3
                elif elmtype == 3: # 4-node quad
                    num_nodes = 4
                elif elmtype == 4: # 4-node tet
                    num_nodes = 4
                else:
                    print ("element type", elmtype, "not implemented")
                    
                nodenums = line[3+numtags:3+numtags+num_nodes]
                nodenums2 = [ pointmap[int(nn)] for nn in nodenums ]


                if elmtype in [2,3]:  # 2d elements

                    # generate facedescriptor for boundary conditon number,
                    # first tag is used as bc-number
                    # element index maps into facedescriptor array
                    if tag1 in facedescriptormap.keys():
                        fdindex = facedescriptormap[tag1]
                    else:
                        fd = FaceDescriptor(bc=tag1)
                        fdindex = mesh.Add(fd)
                        facedescriptormap[tag1] = fdindex
                    
                    mesh.Add (Element2D(fdindex, nodenums2))
                    
                if elmtype == 4:  # 4-node tet
                    nodenums2 = [ pointmap[int(nn)] for nn in nodenums ]
                    mesh.Add (Element3D(tag1, [nodenums2[0],nodenums2[1],nodenums2[3],nodenums2[2]]))
                print (nodenums2)
    # print (mesh)
    # print (pointmap)

    return mesh
    
# ReadGMSH ("test.gmsh")

