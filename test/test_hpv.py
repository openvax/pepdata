from epitopes import hpv

def test_hpv_load_tcell():
    df = hpv.load_tcell()
    assert df is not None

def test_hpv_load_mhc():
    df = hpv.load_mhc()
    assert df is not None