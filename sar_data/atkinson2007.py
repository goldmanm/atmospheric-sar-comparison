#!/usr/bin/env python
# encoding: utf-8

name = "Atkinson2007"
longDesc = u"""
The reaction site *3 needs a lone pair in order to react. It cannot be 2S or 4S.
"""
entry(
    label = "parent",
    group = 
"""
1 *3 R u1 
""",
    data = None
)

entry(
    label = "methyl",
    molecule = """
multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""",
    data = 12.9
)

entry(
    label = "primary_C",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    H u0 {1,S}
4    R!H u0 {1,S}
""",
    data = 9.5
)

entry(
    label = "ethyl",
    molecule = 
"""
multiplicity 2
1 C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 *3 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {1,S}

""",
    data = 9.8
)

entry(
    label = "1_propyl",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {6,S} {7,S} {8,S}
2  C u0 p0 c0 {1,S} {3,S} {9,S} {10,S}
3 *3 C u1 p0 c0 {2,S} {4,S} {5,S}
4  H u0 p0 c0 {3,S}
5  H u0 p0 c0 {3,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {1,S}
9  H u0 p0 c0 {2,S}
10 H u0 p0 c0 {2,S}

""",
    data = 9.6
)

entry(
    label = "primary_alcohol",
    molecule = 
"""
multiplicity 2
1 O u0 p2 c0 {2,S} {5,S}
2 *3 C u1 p0 c0 {1,S} {3,S} {4,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {1,S}

""",
    data = 6.8
)

entry(
    label = "secondary_C",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    H u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
""",
    data = 8.2
)

entry(
    label = "2_propyl",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 *3 C u1 p0 c0 {1,S} {3,S} {4,S}
3  H u0 p0 c0 {2,S}
4  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {4,S}
10 H u0 p0 c0 {4,S}
""",
    data = 8.4
)

entry(
    label = "secondary_alcohol",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    O u0 {1,S} {5,S}
3    R!H u0 {1,S}
4    H u0 {1,S}
5    H u0 {2,S}
""",
    data = 5.2
)

entry(
    label = "tertiary_C",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    R!H u0 {1,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
""",
    data = 7.1
)

entry(
    label = "t_butyl",
    molecule = 
"""
multiplicity 2
1  C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2 *3 C u1 p0 c0 {1,S} {3,S} {4,S}
3  C u0 p0 c0 {2,S} {8,S} {9,S} {10,S}
4  C u0 p0 c0 {2,S} {11,S} {12,S} {13,S}
5  H u0 p0 c0 {1,S}
6  H u0 p0 c0 {1,S}
7  H u0 p0 c0 {1,S}
8  H u0 p0 c0 {3,S}
9  H u0 p0 c0 {3,S}
10 H u0 p0 c0 {3,S}
11 H u0 p0 c0 {4,S}
12 H u0 p0 c0 {4,S}
13 H u0 p0 c0 {4,S}

""",
    data = 7.2
)

entry(
    label = "tertiary_alcohol",
    group = 
"""
1 *3 C u1 {2,S} {3,S} {4,S}
2    O u0 {1,S} {5,S}
3    R!H u0 {1,S}
4    R!H u0 {1,S}
5    H u0 {2,S}
""",
    data = 4.8
)

entry(
    label = "oxy",
    group = 
"""
1 *3 O u1 
""",
    data = 7.5
)

entry(
    label = "formyl",
    molecule = 
"""
multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 O u0 p2 c0 {1,D}
""",
    data = 10.
)

entry(
    label = "carbonyl",
    group = 
"""
1 *3 C u1 {2,D} {3,S}
2    O u0 {1,D}
3    R!H u0 {1,S}
""",
    data = 5.
)

tree(
"""
L1: parent
    L2: methyl
    L2: primary_C
        L3: ethyl
        L3: 1_propyl
        L3: primary_alcohol
    L2: secondary_C
        L3: 2_propyl
        L3: secondary_alcohol
    L2: tertiary_C
        L3: t_butyl
        L3: tertiary_alcohol
    L2: oxy
    L2: formyl
    L2: carbonyl
"""
)

