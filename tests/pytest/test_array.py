from pyngcore import *
from numpy import sort, array

def test_array_numpy():
    a = Array_I_S(5)
    a[:] = 0
    a[3:] = 2
    assert(sum(a) == 4)
    a[1] = 5
    b = sort(a)
    assert(all(b == array([0,0,2,2,5])))
    assert(all(a == array([0,5,0,2,2])))
    a.NumPy().sort()
    assert(all(a == array([0,0,2,2,5])))

if __name__ == "__main__":
    test_array_numpy()
