from netgen.meshing import *


def ReadGmsh(filename):
    meshdim = 1
    with open(filename + ".msh", 'r') as f:
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

    f = open(filename + ".msh", 'r')
    mesh = Mesh(dim=meshdim)

    pointmap = {}
    facedescriptormap = {}
    namemap = {}
    materialmap = {}
    bbcmap = {}

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
                namemap[int(line.split()[1])] = line.split()[2][1:-1]

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

                if elmtype == 1:  # 2-node line
                    num_nodes = 2
                elif elmtype == 2:  # 3-node trig
                    num_nodes = 3
                elif elmtype == 3:  # 4-node quad
                    num_nodes = 4
                elif elmtype == 4:  # 4-node tet
                    num_nodes = 4
                elif elmtype == 5:  # 8-node hex
                    num_nodes = 8
                elif elmtype == 6:  # 6-node prism
                    num_nodes = 6
                elif elmtype == 7:  # 5-node pyramid
                    num_nodes = 5
                elif elmtype == 15:  # 1-node point
                    num_nodes = 1
                else:
                    raise Exception("element type", elmtype, "not implemented")

                nodenums = line[3 + numtags:3 + numtags + num_nodes]
                nodenums2 = [pointmap[int(nn)] for nn in nodenums]

                if elmtype == 1:
                    if meshdim == 3:
                        if tags[1] in bbcmap:
                            index = bbcmap[tags[1]]
                        else:
                            index = len(bbcmap) + 1
                            if len(namemap):
                                mesh.SetCD2Name(index, namemap[tags[0]])
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
                                fd.bcname = namemap[tags[0]]
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
                                mesh.SetMaterial(index, namemap[tags[0]])
                            else:
                                mesh.SetMaterial(index, "line" + str(tags[1]))
                            materialmap[tags[1]] = index

                    mesh.Add(Element1D(index=index, vertices=nodenums2))

                if elmtype in [2, 3]:  # 2d elements
                    if meshdim == 3:
                        if tags[1] in facedescriptormap.keys():
                            index = facedescriptormap[tags[1]]
                        else:
                            index = len(facedescriptormap) + 1
                            fd = FaceDescriptor(bc=index)
                            if len(namemap):
                                fd.bcname = namemap[tags[0]]
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
                                mesh.SetMaterial(index, namemap[tags[0]])
                            else:
                                mesh.SetMaterial(index, "surf" + str(tags[1]))
                            materialmap[tags[1]] = index

                    mesh.Add(Element2D(index, nodenums2))

                if elmtype in [4, 5, 6, 7]:  # volume elements
                    if tags[1] in materialmap:
                        index = materialmap[tags[1]]
                    else:
                        index = len(materialmap) + 1
                        if len(namemap):
                            mesh.SetMaterial(index, namemap[tags[0]])
                        else:
                            mesh.SetMaterial(index, "vol" + str(tags[1]))
                        materialmap[tags[1]] = index

                    nodenums2 = [pointmap[int(nn)] for nn in nodenums]

                    if elmtype == 4:
                        mesh.Add(Element3D(index, [nodenums2[0], nodenums2[1], nodenums2[3], nodenums2[2]]))
                    elif elmtype in [5, 6, 7]:
                        mesh.Add(Element3D(index, nodenums2))

    return mesh
