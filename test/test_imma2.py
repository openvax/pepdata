from epitopes import imma2

def test_imma2_no_overlap():
    """
    IMMA2: Make sure the immunogenic and non-immunogenic sets have no common
    strings
    """
    imm, non = imma2.load_imma2()
    assert len(imm.intersection(non)) == 0

def test_imma2_ngrams_same_size():
    """
    IMMA2: The number of samples in the n-gram dataset should be the same
    as the original sets
    """
    X, Y = imma2.load_imma2_ngrams()
    assert len(X) == len(Y)
    imm, non = imma2.load_imma2()
    assert len(X) == len(imm) + len(non)
