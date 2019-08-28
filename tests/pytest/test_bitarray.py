from pyngcore import BitArray

def test_bitarray():
    a = BitArray(498)
    assert len(a) == 498

    a.Set()
    for b in a:
        assert b == True

    a.Clear(23)
    assert a[22] == True
    assert a[23] == False
    assert a[24] == True

    a.Clear()
    for b in a:
        assert b == False

    a.Set(23)
    assert a[22] == False
    assert a[23] == True
    assert a[24] == False

    a.Clear()
    a[100:200:9] = True
    for i in range(len(a)):
        assert a[i] == bool(100<=i and i<200 and i%9==100%9)

    ac = ~a

    for b,bc in zip(a,ac):
        assert b == (not bc)


