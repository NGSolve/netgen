import os
from sys import path
from sys import platform as __platform

if __platform.startswith('linux'):
    path.append(os.path.dirname(__file__) + '/../../..')
if __platform.startswith('win'):
    path.append(os.path.dirname(__file__) + '/../../../bin')
if __platform.startswith('darwin'):
    path.append(os.path.dirname(__file__) + '/../../../../../MacOS')

import libngpy
del path
