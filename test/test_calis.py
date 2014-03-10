from epitopes import calis

def test_calis_s1():
    """
    Calis S1: Make sure humans are a subset of all species
    """
    df_all = calis.load_s1(human=False)
    df_human = calis.load_s1(human=True)
    assert len(df_all) > len(df_human)

def test_calis_s2():
    """
    Calis S2: Make sure humans are a subset of all species
    """
    df_all = calis.load_s2(human=False)
    df_human = calis.load_s2(human=True)
    assert len(df_all) > len(df_human)

def test_calis_s1_ngrams():
    """
    Calis S1: Number of samples in ngram dataset should be same as strings
    """
    X,Y = calis.load_s1_ngrams()
    imm, non = calis.load_s1_classes()
    assert len(X) == len(Y)
    assert len(X) == len(imm) + len(non)

def test_calis_s2_ngrams():
    """
    Calis S2: Number of samples in ngram dataset should be same as strings
    """
    X,Y = calis.load_s2_ngrams()
    imm, non = calis.load_s2_classes()
    assert len(X) == len(Y)
    assert len(X) == len(imm) + len(non)
