# occupations in format used for the DFT code

occupation_exceptions = {
    # d-block
    24: [[2, 2, 2, 1], [6, 6], [5]],
    29: [[2, 2, 2, 1], [6, 6], [10]],
    41: [[2, 2, 2, 2, 1], [6, 6, 6], [10, 4]],
    42: [[2, 2, 2, 2, 1], [6, 6, 6], [10, 5]],
    44: [[2, 2, 2, 2, 1], [6, 6, 6], [10, 7]],
    45: [[2, 2, 2, 2, 1], [6, 6, 6], [10, 8]],
    46: [[2, 2, 2, 2], [6, 6, 6], [10, 10]],
    47: [[2, 2, 2, 2, 1], [6, 6, 6], [10, 10]],
    78: [[2, 2, 2, 2, 2, 1], [6, 6, 6, 6], [10, 10, 9], [14]],
    79: [[2, 2, 2, 2, 2, 1], [6, 6, 6, 6], [10, 10, 10], [14]],
    103: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6, 1], [10, 10, 10], [14, 14]],
    # f-block
    57: [[2, 2, 2, 2, 2, 2], [6, 6, 6, 6], [10, 10, 1]],
    58: [[2, 2, 2, 2, 2, 2], [6, 6, 6, 6], [10, 10, 1], [1]],
    64: [[2, 2, 2, 2, 2, 2], [6, 6, 6, 6], [10, 10, 1], [7]],
    89: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 1], [14]],
    90: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 2], [14]],
    91: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 1], [14, 2]],
    92: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 1], [14, 3]],
    93: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 1], [14, 4]],
    96: [[2, 2, 2, 2, 2, 2, 2], [6, 6, 6, 6, 6], [10, 10, 10, 1], [14, 7]]
}


def generate_occupations(z):
    assert 0 <= z <= 118

    # exceptions are hard-coded
    if z in occupation_exceptions:
        return occupation_exceptions[z]

    # orbital order per Aufbau principle: (n, l)
    orbital_order = [
        (1, 0),  # 1s
        (2, 0),  # 2s
        (2, 1),  # 2p
        (3, 0),  # 3s
        (3, 1),  # 3p
        (4, 0),  # 4s
        (3, 2),  # 3d
        (4, 1),  # 4p
        (5, 0),  # 5s
        (4, 2),  # 4d
        (5, 1),  # 5p
        (6, 0),  # 6s
        (4, 3),  # 4f
        (5, 2),  # 5d
        (6, 1),  # 6p
        (7, 0),  # 7s
        (5, 3),  # 5f
        (6, 2),  # 6d
        (7, 1),  # 7p
    ]

    # fill orbitals in order until no more electrons are left
    occupations = list()
    for n, l in orbital_order:
        if z == 0:
            break
        if len(occupations) < l + 1:
            occupations.append(list())
        max_e = 2 * (2 * l + 1)
        if z >= max_e:
            occupations[l].append(max_e)
            z -= max_e
        else:
            occupations[l].append(z)
            z = 0
    return occupations


def occupation_string(string):
    occupations = list()
    stack = list(reversed(string.split()))
    while len(stack) > 0:
        s = stack.pop()
        if s[0] == "[":
            stack.extend({"He": ["1s2"],
                          "Ne": ["2p6", "2s2", "[He]"],
                          "Ar": ["3p6", "3s2", "[Ne]"],
                          "Kr": ["4p6", "3d10", "4s2", "[Ar]"],
                          "Xe": ["5p6", "4d10", "5s2", "[Kr]"],
                          "Rn": ["6p6", "5d10", "4f14", "6s2", "[Xe]"]}[s[1:3]])
        else:
            n = int(s[0])
            l = {"s": 0, "p": 1, "d": 2, "f": 3}[s[1]]
            e = int(s[2:])
            while len(occupations) < l + 1:
                occupations.append(list())
            while len(occupations[l]) < n - l:
                occupations[l].append(0)
            occupations[l][n - l - 1] = e
    return occupations


# copied from ASE (ase.data.__init__.py)

chemical_symbols = [
    # 0
    'X',
    # 1
    'H', 'He',
    # 2
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    # 3
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    # 4
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    # 5
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    # 6
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    # 7
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']

atomic_numbers = {symbol: Z for Z, symbol in enumerate(chemical_symbols)}
