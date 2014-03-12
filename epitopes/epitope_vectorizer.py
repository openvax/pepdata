from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

from reduced_alphabet import make_alphabet_transformer

def make_count_vectorizer(reduced_alphabet, max_ngram):
    if reduced_alphabet is None:
        preprocessor = None
    else:
        preprocessor = make_alphabet_transformer(reduced_alphabet)

    return CountVectorizer(
        analyzer='char',
        ngram_range=(1, max_ngram),
        dtype=np.float,
        preprocessor = preprocessor)

class EpitopeVectorizer(object):
    def __init__(
            self,
            max_ngram = 1,
            normalize_row = True,
            reduced_alphabet = None,
            training_already_reduced = False):
        self.reduced_alphabet = reduced_alphabet
        self.max_ngram = max_ngram
        self.normalize_row = normalize_row
        self.training_already_reduced = training_already_reduced

    def fit_transform(self, amino_acid_strings):
        self.count_vectorizer = \
            make_count_vectorizer(self.reduced_alphabet, self.max_ngram)

        if self.training_already_reduced:
            c = make_count_vectorizer(None, self.max_ngram)
        else:
            c = self.count_vectorizer

        X = c.fit_transform(amino_acid_strings).todense()
        if self.normalize_row:
            X = normalize(X, norm='l1')

    def fit(self, amino_acid_strings):
        self.fit_transform(amino_acid_strings)

    def transform(self, amino_acid_strings):
        X = self.count_vectorizer.transform(amino_acid_strings).todense()
        if self.normalize_row:
            X = normalize(X, norm='l1')
        return X
