# Netgen mesh generator

NETGEN is an automatic 3d tetrahedral mesh generator. It accepts input from constructive solid geometry (CSG) or boundary representation (BRep) from STL file format. The connection to a geometry kernel allows the handling of IGES and STEP files. NETGEN contains modules for mesh optimization and hierarchical mesh refinement. Netgen 6.x supports scripting via a Python interface. Netgen is open source based on the LGPL license. It is available for Unix/Linux, Windows, and OSX.

Find the Open Source Community on https://ngsolve.org
Support & Services: https://cerbsim.com

## Installing Locally

After installing all dependencies, run `git clone <netgen-repo> && pip install ./netgen`

### System Dependencies

#### Ubuntu

```
sudo apt install python3 git cmake tcl-dev tk-dev libglu1-mesa-dev libxmu-dev
```

### Python Dependencies

```
pip install scikit-build requests netgen-occt-devel
```
