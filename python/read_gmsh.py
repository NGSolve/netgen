from netgen.meshing import *

def ReadGmsh(filename):
    if not filename.endswith(".msh"):
        filename += ".msh"
    meshdim = 1
    with open(filename, 'r') as f:
        while f.readline().split()[0] != "$Elements":
            pass
        nelem = int(f.readline())
        for i in range(nelem):
            line = f.readline()
            eltype = int(line.split()[1])
            if eltype > 1 and eltype != 15:
                meshdim = 2
            if eltype > 3 and eltype != 15:
                meshdim = 3
                break

    f = open(filename, 'r')
    mesh = Mesh(dim=meshdim)

    pointmap = {}
    facedescriptormap = {}
    namemap = { 0 : { 0 : "default" },
                1: { 0 : "default" },
                2: { 0 : "default" },
                3: { 0 : "default" } }
    materialmap = {}
    bbcmap = {}

    segm = 1
    trig = 2
    quad = 3
    tet = 4
    hex = 5
    prism = 6
    pyramid = 7
    segm3 = 8      # 2nd order line
    trig6 = 9      # 2nd order trig
    tet10 = 11     # 2nd order tet
    point = 15
    quad8 = 16     # 2nd order quad
    hex20 = 17     # 2nd order hex
    prism15 = 18   # 2nd order prism
    pyramid13 = 19 # 2nd order pyramid
    segms = [segm, segm3]
    trigs = [trig, trig6]
    quads = [quad, quad8]
    tets = [tet, tet10]
    hexes = [hex, hex20]
    prisms = [prism, prism15]
    pyramids = [pyramid, pyramid13]
    elem0d = [point]
    elem1d = segms
    elem2d = trigs + quads
    elem3d = tets + hexes + prisms + pyramids

    num_nodes_map = { segm : 2,
                      trig : 3,
                      quad : 4,
                      tet : 4,
                      hex : 8,
                      prism : 6,
                      pyramid : 5,
                      segm3 : 3,
                      trig6 : 6,
                      tet10 : 10,
                      point : 1,
                      quad8 : 8,
                      hex20 : 20,
                      prism15 : 18,
                      pyramid13 : 19 }

    while True:
        line = f.readline()
        if line == "":
            break

        if line.split()[0] == "$PhysicalNames":
            print('WARNING: Physical groups detected - Be sure to define them for every geometrical entity.')
            numnames = int(f.readline())
            for i in range(numnames):
                f.readline
                line = f.readline()
                namemap[int(line.split()[0])][int(line.split()[1])] = line.split()[2][1:-1]

        if line.split()[0] == "$Nodes":
            num = int(f.readline().split()[0])
            for i in range(num):
                line = f.readline()
                nodenum, x, y, z = line.split()[0:4]
                pnum = mesh.Add(MeshPoint(Pnt(float(x), float(y), float(z))))
                pointmap[int(nodenum)] = pnum

        if line.split()[0] == "$Elements":
            num = int(f.readline().split()[0])

            for i in range(num):
                line = f.readline().split()
                elmnum = int(line[0])
                elmtype = int(line[1])
                numtags = int(line[2])
                # the first tag is the physical group nr, the second tag is the group nr of the dim
                tags = [int(line[3 + k]) for k in range(numtags)]

                if elmtype not in num_nodes_map:
                    raise Exception("element type", elmtype, "not implemented")
                num_nodes = num_nodes_map[elmtype]

                nodenums = line[3 + numtags:3 + numtags + num_nodes]
                nodenums2 = [pointmap[int(nn)] for nn in nodenums]

                if elmtype in elem1d:
                    if meshdim == 3:
                        if tags[1] in bbcmap:
                            index = bbcmap[tags[1]]
                        else:
                            index = len(bbcmap) + 1
                            if len(namemap):
                                mesh.SetCD2Name(index, namemap[1][tags[0]])
                            else:
                                mesh.SetCD2Name(index, "line" + str(tags[1]))
                            bbcmap[tags[1]] = index

                    elif meshdim == 2:
                        if tags[1] in facedescriptormap.keys():
                            index = facedescriptormap[tags[1]]
                        else:
                            index = len(facedescriptormap) + 1
                            fd = FaceDescriptor(bc=index)
                            if len(namemap):
                                fd.bcname = namemap[1][tags[0]]
                            else:
                                fd.bcname = 'line' + str(tags[1])
                            mesh.SetBCName(index - 1, fd.bcname)
                            mesh.Add(fd)
                            facedescriptormap[tags[1]] = index
                    else:
                        if tags[1] in materialmap:
                            index = materialmap[tags[1]]
                        else:
                            index = len(materialmap) + 1
                            if len(namemap):
                                mesh.SetMaterial(index, namemap[1][tags[0]])
                            else:
                                mesh.SetMaterial(index, "line" + str(tags[1]))
                            materialmap[tags[1]] = index

                    mesh.Add(Element1D(index=index, vertices=nodenums2))

                if elmtype in elem2d:  # 2d elements
                    if meshdim == 3:
                        if tags[1] in facedescriptormap.keys():
                            index = facedescriptormap[tags[1]]
                        else:
                            index = len(facedescriptormap) + 1
                            fd = FaceDescriptor(bc=index)
                            if len(namemap):
                                fd.bcname = namemap[2][tags[0]]
                            else:
                                fd.bcname = "surf" + str(tags[1])
                            mesh.SetBCName(index - 1, fd.bcname)
                            mesh.Add(fd)
                            facedescriptormap[tags[1]] = index
                    else:
                        if tags[1] in materialmap:
                            index = materialmap[tags[1]]
                        else:
                            index = len(materialmap) + 1
                            if len(namemap):
                                mesh.SetMaterial(index, namemap[2][tags[0]])
                            else:
                                mesh.SetMaterial(index, "surf" + str(tags[1]))
                            materialmap[tags[1]] = index

                    if elmtype in trigs:
                        ordering = [i for i in range(3)]
                        if elmtype == trig6:
                            ordering += [4,5,3]
                    if elmtype in quads:
                        ordering = [i for i in range(4)]
                        if elmtype == quad8:
                            ordering += [4, 6, 7, 5]
                    mesh.Add(Element2D(index, [nodenums2[i] for i in ordering]))

                if elmtype in elem3d:  # volume elements
                    if tags[1] in materialmap:
                        index = materialmap[tags[1]]
                    else:
                        index = len(materialmap) + 1
                        if len(namemap):
                            mesh.SetMaterial(index, namemap[3][tags[0]])
                        else:
                            mesh.SetMaterial(index, "vol" + str(tags[1]))
                        materialmap[tags[1]] = index

                    nodenums2 = [pointmap[int(nn)] for nn in nodenums]

                    if elmtype in tets:
                        ordering = [0,1,2,3]
                        if elmtype == tet10:
                            ordering += [4,6,7,5,9,8]
                    elif elmtype in hexes:
                        ordering = [0,1,5,4,3,2,6,7]
                        if elmtype == hex20:
                            ordering += [8,16,10,12,13,19,15,14,9,11,18,17]
                    elif elmtype in prisms:
                        ordering = [0,2,1,3,5,4]
                        if elmtype == prism15:
                            ordering += [7,6,9,8,11,10,13,12,14]
                    elif elmtype in pyramids:
                        ordering = [3,2,1,0,4]
                        if elmtype == pyramid13:
                            ordering += [10,5,6,8,12,11,9,7]
                    mesh.Add(Element3D(index, [nodenums2[i] for i in ordering]))

    return mesh
