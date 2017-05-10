import os
from sys import path
from sys import platform as __platform

if __platform.startswith('linux'):
    _netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'../../../../bin'))
    _netgen_lib_dir=os.path.realpath(os.path.join(os.path.dirname(__file__),'../../../../lib'))
if __platform.startswith('win'):
    _netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__) + '/../../../bin'))
    _netgen_lib_dir=_netgen_bin_dir
if __platform.startswith('darwin'):
    _netgen_bin_dir=os.path.realpath(os.path.join(os.path.dirname(__file__) + '/../../../../../MacOS'))
    _netgen_lib_dir=_netgen_bin_dir

path.append(_netgen_lib_dir)
path.append(_netgen_bin_dir)

import libngpy
del path
