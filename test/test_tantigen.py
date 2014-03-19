from epitopes import tantigen

def test_load_tcell():
    df = tantigen.load_tcell()
    assert df is not None

def test_load_mhc():
    df = tantigen.load_mhc()
    assert df is not None
