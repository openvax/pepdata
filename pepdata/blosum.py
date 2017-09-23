
def parse_blosum_table(table, coeff_type=int, key_type='row'):
    """
    Parse a table of pairwise amino acid coefficient (e.g. BLOSUM50)
    """

    lines = table.split("\n")
    # drop comments
    lines = [line for line in lines if not line.startswith("#")]
    # drop CR endline characters
    lines = [line.replace("\r", "") for line in lines]
    # skip empty lines
    lines = [line for line in lines if line]

    labels = lines[0].split()

    if len(labels) < 20:
        raise ValueError(
            "Expected 20+ amino acids but first line '%s' has %d fields" % (
                lines[0],
                len(labels)))
    coeffs = {}
    for line in lines[1:]:

        fields = line.split()
        assert len(fields) >= 21, \
            "Expected AA and 20+ coefficients but '%s' has %d fields" % (
                line, len(fields))
        x = fields[0]
        for i, coeff_str in enumerate(fields[1:]):
            y = labels[i]
            coeff = coeff_type(coeff_str)
            if key_type == 'pair':
                coeffs[(x, y)] = coeff
            elif key_type == 'pair_string':
                coeffs[x + y] = coeff
            else:
                assert key_type == 'row', "Unknown key type: %s" % key_type
                if x not in coeffs:
                    coeffs[x] = {}
                coeffs[x][y] = coeff
    return coeffs


with open(join(MATRIX_DIR, 'BLOSUM30'), 'r') as f:
    blosum30 = parse_blosum_table(f.read())

with open(join(MATRIX_DIR, 'BLOSUM50'), 'r') as f:
    blosum50 = parse_blosum_table(f.read())

with open(join(MATRIX_DIR, 'BLOSUM62'), 'r') as f:
    blosum62 = parse_blosum_table(f.read())
