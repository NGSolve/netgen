from .pyngcore import *

# <size_t> is the same as <unsigned int> on 32 bit arches
# in which case Array_I_S is not defined by python_ngcore_export.cpp.
# In this case identify it with Array_I_U.
try: Array_I_S
except NameError: Array_I_S=Array_I_U
