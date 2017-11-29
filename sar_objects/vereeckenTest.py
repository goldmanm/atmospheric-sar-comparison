# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:47:20 2017

@author: mark
"""
import unittest

from vereecken2009 import sar, get_atom_groups, vereecken_get_groups, get_cyclic_groups

from rmgpy.thermo.nasa import NASA, NASAPolynomial
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.data.kinetics.family import TemplateReaction, KineticsFamily
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.thermo import ThermoData

class SarRateEstimationTest(unittest.TestCase):

    def setUp(self):
        self.sar = sar
        # ethoxy reaction
        self.reaction = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([35.4385,39.2459,43.7646,48.1997,55.9401,61.965,61.965],'J/(mol*K)'), H298=(-108.575,'kJ/mol'), S298=(218.834,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH)"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,S} {4,D}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4 *2 O u0 p2 c0 {1,D}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([9.331,9.988,10.611,11.274,12.457,13.427,15.105],'cal/(mol*K)'), H298=(34.893,'kcal/mol'), S298=(46.3704,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([61.879,72.3026,83.2296,93.9954,112.876,127.558,150.145],'J/(mol*K)'), H298=(-17.2758,'kJ/mol'), S298=(284.871,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 *1 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3 *2 O u1 p2 c0 {2,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {2,S}
""")])], pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([35.4385,39.2459,43.7646,48.1997,55.9401,61.965,61.965],'J/(mol*K)'), H298=(-108.575,'kJ/mol'), S298=(218.834,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH)"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,S} {4,D}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4 *2 O u0 p2 c0 {1,D}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([61.879,72.3026,83.2296,93.9954,112.876,127.558,150.145],'J/(mol*K)'), H298=(-17.2758,'kJ/mol'), S298=(284.871,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 *1 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3 *2 O u1 p2 c0 {2,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {2,S}
""")])], [Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([9.331,9.988,10.611,11.274,12.457,13.427,15.105],'cal/(mol*K)'), H298=(34.893,'kcal/mol'), S298=(46.3704,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([61.879,72.3026,83.2296,93.9954,112.876,127.558,150.145],'J/(mol*K)'), H298=(-17.2758,'kJ/mol'), S298=(284.871,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 *1 C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3 *2 O u1 p2 c0 {2,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {2,S}
""")])]], family='R_Addition_MultipleBond', template=['Cd_R', 'CsJ'])
        # ipropyloxy
        self.degenerate_reaction = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([51.55,63.29,74.38,84.19,100.1,112.02,130.2],'J/(mol*K)'), H298=(-166.3,'kJ/mol'), S298=(263.446,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Od)HHH) + group(Cds-OdCsH)"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *1 C u0 p0 c0 {1,S} {6,D} {7,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6 *2 O u0 p2 c0 {2,D}
7    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([9.331,9.988,10.611,11.274,12.457,13.427,15.105],'cal/(mol*K)'), H298=(34.893,'kcal/mol'), S298=(46.3704,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([87.5579,108.135,126.429,141.841,165.681,183.615,211.679],'J/(mol*K)'), H298=(-46.0526,'kJ/mol'), S298=(302.672,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  *3 C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  *2 O u1 p2 c0 {1,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
""")])], degeneracy=2.0, pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([51.55,63.29,74.38,84.19,100.1,112.02,130.2],'J/(mol*K)'), H298=(-166.3,'kJ/mol'), S298=(263.446,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Od)HHH) + group(Cds-OdCsH)"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *1 C u0 p0 c0 {1,S} {6,D} {7,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6 *2 O u0 p2 c0 {2,D}
7    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([87.5579,108.135,126.429,141.841,165.681,183.615,211.679],'J/(mol*K)'), H298=(-46.0526,'kJ/mol'), S298=(302.672,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  *3 C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  *2 O u1 p2 c0 {1,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
""")])], [Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([9.331,9.988,10.611,11.274,12.457,13.427,15.105],'cal/(mol*K)'), H298=(34.893,'kcal/mol'), S298=(46.3704,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2    H u0 p0 c0 {1,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([87.5579,108.135,126.429,141.841,165.681,183.615,211.679],'J/(mol*K)'), H298=(-46.0526,'kJ/mol'), S298=(302.672,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3  *3 C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  *2 O u1 p2 c0 {1,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
""")])]], family='R_Addition_MultipleBond', template=['Cd_R', 'CsJ'])
        # C=CC=O + [O]C(C)C
        self.multifunc_rxn = TemplateReaction(reactants=[Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.27135,0.0262311,-9.29123e-06,-4.78373e-09,3.34805e-12,-9335.73,19.4981], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.81119,0.0171143,-7.48342e-06,1.42522e-09,-9.17468e-14,-10784.1,-4.8588], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(298,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""C2H3CHO""", comment="""Thermo library: JetSurF2.0"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,D} {3,S} {4,S}
2    C u0 p0 c0 {1,D} {5,S} {6,S}
3 *1 C u0 p0 c0 {1,S} {7,S} {8,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8 *2 O u0 p2 c0 {3,D}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([20.827,25.556,29.742,33.303,38.965,43.242,50],'cal/(mol*K)','+|-',[1.4,1.6,1.7,1.8,1.8,1.7,1.4]), H298=(-10.664,'kcal/mol','+|-',0.9), S298=(72.525,'cal/(mol*K)','+|-',1.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""ipropoxy""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  *3 O u1 p2 c0 {1,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([160.168,192.997,225.501,252.298,294.769,327.347,372.12],'J/(mol*K)'), H298=(-152.072,'kJ/mol'), S298=(456.011,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  *1 C u0 p0 c0 {5,S} {7,S} {9,S} {16,S}
3     C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4     C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5     C u0 p0 c0 {2,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7  *3 O u0 p2 c0 {1,S} {2,S}
8     H u0 p0 c0 {1,S}
9  *2 O u1 p2 c0 {2,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {2,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])], pairs=[[Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.27135,0.0262311,-9.29123e-06,-4.78373e-09,3.34805e-12,-9335.73,19.4981], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.81119,0.0171143,-7.48342e-06,1.42522e-09,-9.17468e-14,-10784.1,-4.8588], Tmin=(1000,'K'), Tmax=(5000,'K'))], Tmin=(298,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""C2H3CHO""", comment="""Thermo library: JetSurF2.0"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,D} {3,S} {4,S}
2    C u0 p0 c0 {1,D} {5,S} {6,S}
3 *1 C u0 p0 c0 {1,S} {7,S} {8,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8 *2 O u0 p2 c0 {3,D}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([160.168,192.997,225.501,252.298,294.769,327.347,372.12],'J/(mol*K)'), H298=(-152.072,'kJ/mol'), S298=(456.011,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  *1 C u0 p0 c0 {5,S} {7,S} {9,S} {16,S}
3     C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4     C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5     C u0 p0 c0 {2,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7  *3 O u0 p2 c0 {1,S} {2,S}
8     H u0 p0 c0 {1,S}
9  *2 O u1 p2 c0 {2,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {2,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])], [Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([20.827,25.556,29.742,33.303,38.965,43.242,50],'cal/(mol*K)','+|-',[1.4,1.6,1.7,1.8,1.8,1.7,1.4]), H298=(-10.664,'kcal/mol','+|-',0.9), S298=(72.525,'cal/(mol*K)','+|-',1.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""ipropoxy""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4  *3 O u1 p2 c0 {1,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([160.168,192.997,225.501,252.298,294.769,327.347,372.12],'J/(mol*K)'), H298=(-152.072,'kJ/mol'), S298=(456.011,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + longDistanceInteraction_noncyclic(OsCs-ST) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {3,S} {4,S} {7,S} {8,S}
2  *1 C u0 p0 c0 {5,S} {7,S} {9,S} {16,S}
3     C u0 p0 c0 {1,S} {10,S} {11,S} {12,S}
4     C u0 p0 c0 {1,S} {13,S} {14,S} {15,S}
5     C u0 p0 c0 {2,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7  *3 O u0 p2 c0 {1,S} {2,S}
8     H u0 p0 c0 {1,S}
9  *2 O u1 p2 c0 {2,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {2,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])]], family='R_Addition_MultipleBond', template=['CO-CdH_O', 'OJ-Cs'])
        # C=C[CH2] and O=COCC
        self.ether_rxn = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([21.312,25.646,29.989,34,40.56,45.291,52.139],'cal/(mol*K)'), H298=(-92.804,'kcal/mol'), S298=(78.538,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""C3H6O2b""", comment="""Thermo library: CHO"""), molecule=[Molecule().fromAdjacencyList("""1     C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2     C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  *1 C u0 p0 c0 {4,S} {10,D} {11,S}
4     O u0 p2 c0 {1,S} {3,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {2,S}
10 *2 O u0 p2 c0 {3,D}
11    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.36318,0.0198138,1.24971e-05,-3.33556e-08,1.58466e-11,19245.6,17.1732], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.50079,0.0143247,-5.67816e-06,1.10808e-09,-9.03639e-14,17482.4,-11.2431], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""aC3H5""", comment="""Thermo library: JetSurF2.0"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,D} {3,S} {4,S}
2    C u0 p0 c0 {1,D} {5,S} {6,S}
3 *3 C u1 p0 c0 {1,S} {7,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8    H u0 p0 c0 {3,S}
"""), Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,S} {3,D} {4,S}
2    C u1 p0 c0 {1,S} {5,S} {6,S}
3 *3 C u0 p0 c0 {1,D} {7,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8    H u0 p0 c0 {3,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([158.711,191.767,224.212,250.959,293.808,326.806,372.067],'J/(mol*K)'), H298=(-140.603,'kJ/mol'), S298=(468.613,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  *1 C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
3     C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
4     C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5     C u0 p0 c0 {1,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7     O u0 p2 c0 {2,S} {3,S}
8     H u0 p0 c0 {1,S}
9     H u0 p0 c0 {1,S}
10 *2 O u1 p2 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {3,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {4,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])], pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([21.312,25.646,29.989,34,40.56,45.291,52.139],'cal/(mol*K)'), H298=(-92.804,'kcal/mol'), S298=(78.538,'cal/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), label="""C3H6O2b""", comment="""Thermo library: CHO"""), molecule=[Molecule().fromAdjacencyList("""1     C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2     C u0 p0 c0 {1,S} {7,S} {8,S} {9,S}
3  *1 C u0 p0 c0 {4,S} {10,D} {11,S}
4     O u0 p2 c0 {1,S} {3,S}
5     H u0 p0 c0 {1,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {2,S}
10 *2 O u0 p2 c0 {3,D}
11    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([158.711,191.767,224.212,250.959,293.808,326.806,372.067],'J/(mol*K)'), H298=(-140.603,'kJ/mol'), S298=(468.613,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  *1 C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
3     C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
4     C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5     C u0 p0 c0 {1,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7     O u0 p2 c0 {2,S} {3,S}
8     H u0 p0 c0 {1,S}
9     H u0 p0 c0 {1,S}
10 *2 O u1 p2 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {3,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {4,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])], [Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[1.36318,0.0198138,1.24971e-05,-3.33556e-08,1.58466e-11,19245.6,17.1732], Tmin=(298,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[6.50079,0.0143247,-5.67816e-06,1.10808e-09,-9.03639e-14,17482.4,-11.2431], Tmin=(1000,'K'), Tmax=(3000,'K'))], Tmin=(298,'K'), Tmax=(3000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""aC3H5""", comment="""Thermo library: JetSurF2.0"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,D} {3,S} {4,S}
2    C u0 p0 c0 {1,D} {5,S} {6,S}
3 *3 C u1 p0 c0 {1,S} {7,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8    H u0 p0 c0 {3,S}
"""), Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,S} {3,D} {4,S}
2    C u1 p0 c0 {1,S} {5,S} {6,S}
3 *3 C u0 p0 c0 {1,D} {7,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {3,S}
8    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([158.711,191.767,224.212,250.959,293.808,326.806,372.067],'J/(mol*K)'), H298=(-140.603,'kJ/mol'), S298=(468.613,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Os-CsCs) + group(Os-CsH) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {5,S} {8,S} {9,S}
2  *1 C u0 p0 c0 {1,S} {7,S} {10,S} {11,S}
3     C u0 p0 c0 {4,S} {7,S} {12,S} {13,S}
4     C u0 p0 c0 {3,S} {14,S} {15,S} {16,S}
5     C u0 p0 c0 {1,S} {6,D} {17,S}
6     C u0 p0 c0 {5,D} {18,S} {19,S}
7     O u0 p2 c0 {2,S} {3,S}
8     H u0 p0 c0 {1,S}
9     H u0 p0 c0 {1,S}
10 *2 O u1 p2 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {3,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {4,S}
16    H u0 p0 c0 {4,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
""")])]], family='R_Addition_MultipleBond', template=['CO-NdH_O', 'CsJ-CdHH'])
        self.gamma_alcohol_beta_alkyl_rxn = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([17.667,21.961,25.631,28.398,31.96,34.133,37.211],'cal/(mol*K)','+|-',[1.6,1.9,1.7,1.3,0.8,0.5,0.3]), H298=(-76.04,'kcal/mol','+|-',0.9), S298=(68.745,'cal/(mol*K)','+|-',1.2), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""HOCH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *1 C u0 p0 c0 {1,S} {6,D} {7,S}
3    O u0 p2 c0 {1,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6 *2 O u0 p2 c0 {2,D}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[9.68072,0.00707718,4.53302e-05,-3.61952e-08,8.06654e-12,1051.17,-19.4358], Tmin=(200,'K'), Tmax=(1375,'K')), NASAPolynomial(coeffs=[11.7983,0.035958,-1.44109e-05,2.46308e-09,-1.52118e-13,-3861.01,-44.1297], Tmin=(1375,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), label="""NEOC5H11""", comment="""Thermo library: CurranPentane"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4     C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  *3 C u1 p0 c0 {1,S} {15,S} {16,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {5,S}
16    H u0 p0 c0 {5,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([198.833,247.364,289.371,324.031,376.268,413.955,471.02],'J/(mol*K)'), H298=(-300.07,'kJ/mol'), S298=(472.927,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  *3 C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  *1 C u0 p0 c0 {2,S} {4,S} {9,S} {12,S}
4     C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5     C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6     C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
7     C u0 p0 c0 {1,S} {21,S} {22,S} {23,S}
8     O u0 p2 c0 {4,S} {24,S}
9  *2 O u1 p2 c0 {3,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {5,S}
16    H u0 p0 c0 {5,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {6,S}
21    H u0 p0 c0 {7,S}
22    H u0 p0 c0 {7,S}
23    H u0 p0 c0 {7,S}
24    H u0 p0 c0 {8,S}
""")])], pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([17.667,21.961,25.631,28.398,31.96,34.133,37.211],'cal/(mol*K)','+|-',[1.6,1.9,1.7,1.3,0.8,0.5,0.3]), H298=(-76.04,'kcal/mol','+|-',0.9), S298=(68.745,'cal/(mol*K)','+|-',1.2), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), label="""HOCH2CHO""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *1 C u0 p0 c0 {1,S} {6,D} {7,S}
3    O u0 p2 c0 {1,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6 *2 O u0 p2 c0 {2,D}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([198.833,247.364,289.371,324.031,376.268,413.955,471.02],'J/(mol*K)'), H298=(-300.07,'kJ/mol'), S298=(472.927,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  *3 C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  *1 C u0 p0 c0 {2,S} {4,S} {9,S} {12,S}
4     C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5     C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6     C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
7     C u0 p0 c0 {1,S} {21,S} {22,S} {23,S}
8     O u0 p2 c0 {4,S} {24,S}
9  *2 O u1 p2 c0 {3,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {5,S}
16    H u0 p0 c0 {5,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {6,S}
21    H u0 p0 c0 {7,S}
22    H u0 p0 c0 {7,S}
23    H u0 p0 c0 {7,S}
24    H u0 p0 c0 {8,S}
""")])], [Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[9.68072,0.00707718,4.53302e-05,-3.61952e-08,8.06654e-12,1051.17,-19.4358], Tmin=(200,'K'), Tmax=(1375,'K')), NASAPolynomial(coeffs=[11.7983,0.035958,-1.44109e-05,2.46308e-09,-1.52118e-13,-3861.01,-44.1297], Tmin=(1375,'K'), Tmax=(5000,'K'))], Tmin=(200,'K'), Tmax=(5000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), label="""NEOC5H11""", comment="""Thermo library: CurranPentane"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2     C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {9,S} {10,S} {11,S}
4     C u0 p0 c0 {1,S} {12,S} {13,S} {14,S}
5  *3 C u1 p0 c0 {1,S} {15,S} {16,S}
6     H u0 p0 c0 {2,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11    H u0 p0 c0 {3,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {5,S}
16    H u0 p0 c0 {5,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([198.833,247.364,289.371,324.031,376.268,413.955,471.02],'J/(mol*K)'), H298=(-300.07,'kJ/mol'), S298=(472.927,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(552.912,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-SQ) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Os-CsH) + group(Os-CsH) + radical(CC(C)OJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1     C u0 p0 c0 {2,S} {5,S} {6,S} {7,S}
2  *3 C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}
3  *1 C u0 p0 c0 {2,S} {4,S} {9,S} {12,S}
4     C u0 p0 c0 {3,S} {8,S} {13,S} {14,S}
5     C u0 p0 c0 {1,S} {15,S} {16,S} {17,S}
6     C u0 p0 c0 {1,S} {18,S} {19,S} {20,S}
7     C u0 p0 c0 {1,S} {21,S} {22,S} {23,S}
8     O u0 p2 c0 {4,S} {24,S}
9  *2 O u1 p2 c0 {3,S}
10    H u0 p0 c0 {2,S}
11    H u0 p0 c0 {2,S}
12    H u0 p0 c0 {3,S}
13    H u0 p0 c0 {4,S}
14    H u0 p0 c0 {4,S}
15    H u0 p0 c0 {5,S}
16    H u0 p0 c0 {5,S}
17    H u0 p0 c0 {5,S}
18    H u0 p0 c0 {6,S}
19    H u0 p0 c0 {6,S}
20    H u0 p0 c0 {6,S}
21    H u0 p0 c0 {7,S}
22    H u0 p0 c0 {7,S}
23    H u0 p0 c0 {7,S}
24    H u0 p0 c0 {8,S}
""")])]], family='R_Addition_MultipleBond', template=['CO-CsH_O', 'CsJ-CsHH'])

        self.formic_acid_reaction = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([11.017,13.147,15.139,16.849,19.423,21.14,23.318],'cal/(mol*K)','+|-',[0.7,0.9,1.1,1.1,0.9,0.7,0.2]), H298=(-90.466,'kcal/mol','+|-',0.1), S298=(59.464,'cal/(mol*K)','+|-',0.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,D} {4,S}
2    O u0 p2 c0 {1,S} {5,S}
3 *2 O u0 p2 c0 {1,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    O u0 p2 c0 {2,S} {3,S}
2 *3 C u1 p0 c0 {1,S} {4,D}
3    H u0 p0 c0 {1,S}
4    O u0 p2 c0 {2,D}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([93.246,108.196,122.511,132.447,148.077,161.173,176.393],'J/(mol*K)'), H298=(-541.543,'kJ/mol'), S298=(377.342,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cds-OdCsOs) + group(Os-CsH) + group(Os-CsH) + group(Os-(Cds-Od)H) + radical(C=OCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 *3 C u0 p0 c0 {1,S} {4,S} {7,D}
3    O u0 p2 c0 {1,S} {8,S}
4    O u0 p2 c0 {2,S} {9,S}
5 *2 O u1 p2 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    O u0 p2 c0 {2,D}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {4,S}
""")])], pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([11.017,13.147,15.139,16.849,19.423,21.14,23.318],'cal/(mol*K)','+|-',[0.7,0.9,1.1,1.1,0.9,0.7,0.2]), H298=(-90.466,'kcal/mol','+|-',0.1), S298=(59.464,'cal/(mol*K)','+|-',0.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,D} {4,S}
2    O u0 p2 c0 {1,S} {5,S}
3 *2 O u0 p2 c0 {1,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([93.246,108.196,122.511,132.447,148.077,161.173,176.393],'J/(mol*K)'), H298=(-541.543,'kJ/mol'), S298=(377.342,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cds-OdCsOs) + group(Os-CsH) + group(Os-CsH) + group(Os-(Cds-Od)H) + radical(C=OCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 *3 C u0 p0 c0 {1,S} {4,S} {7,D}
3    O u0 p2 c0 {1,S} {8,S}
4    O u0 p2 c0 {2,S} {9,S}
5 *2 O u1 p2 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    O u0 p2 c0 {2,D}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {4,S}
""")])], [Species(label="", thermo=NASA(polynomials=[NASAPolynomial(coeffs=[2.92208,0.00762454,3.29884e-06,-1.07135e-08,5.11587e-12,-23028.2,11.2926], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[5.39206,0.00411221,-1.48195e-06,2.39875e-10,-1.43903e-14,-23860.7,-2.23529], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""HOCO""", comment="""Thermo library: FFCM1(-)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    O u0 p2 c0 {2,S} {3,S}
2 *3 C u1 p0 c0 {1,S} {4,D}
3    H u0 p0 c0 {1,S}
4    O u0 p2 c0 {2,D}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([93.246,108.196,122.511,132.447,148.077,161.173,176.393],'J/(mol*K)'), H298=(-541.543,'kJ/mol'), S298=(377.342,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsOsOsH) + group(Cds-OdCsOs) + group(Os-CsH) + group(Os-CsH) + group(Os-(Cds-Od)H) + radical(C=OCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1 *1 C u0 p0 c0 {2,S} {3,S} {5,S} {6,S}
2 *3 C u0 p0 c0 {1,S} {4,S} {7,D}
3    O u0 p2 c0 {1,S} {8,S}
4    O u0 p2 c0 {2,S} {9,S}
5 *2 O u1 p2 c0 {1,S}
6    H u0 p0 c0 {1,S}
7    O u0 p2 c0 {2,D}
8    H u0 p0 c0 {3,S}
9    H u0 p0 c0 {4,S}
""")])]], family='R_Addition_MultipleBond', template=['CO-NdH_O', 'CO_rad/NonDe'])
        self.beta_propyl_cyclic_reaction = TemplateReaction(reactants=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([11.017,13.147,15.139,16.849,19.423,21.14,23.318],'cal/(mol*K)','+|-',[0.7,0.9,1.1,1.1,0.9,0.7,0.2]), H298=(-90.466,'kcal/mol','+|-',0.1), S298=(59.464,'cal/(mol*K)','+|-',0.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,D} {4,S}
2    O u0 p2 c0 {1,S} {5,S}
3 *2 O u0 p2 c0 {1,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([13.614,17.846,21.434,24.32,28.647,31.786,36.685],'cal/(mol*K)','+|-',[1.3,1.6,1.6,1.6,1.5,1.4,1]), H298=(69.882,'kcal/mol','+|-',0.9), S298=(61.584,'cal/(mol*K)','+|-',0.7), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""cC3H5""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2    C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3 *3 C u1 p0 c0 {1,S} {2,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
""")])], products=[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([103.524,128.448,152.586,171.087,199.997,222.111,252.569],'J/(mol*K)'), H298=(-106.461,'kJ/mol'), S298=(371.953,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Os-CsH) + group(Os-CsH) + ring(Cyclopropane) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2     C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  *1 C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5     O u0 p2 c0 {4,S} {13,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11 *2 O u1 p2 c0 {4,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
""")])], pairs=[[Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([11.017,13.147,15.139,16.849,19.423,21.14,23.318],'cal/(mol*K)','+|-',[0.7,0.9,1.1,1.1,0.9,0.7,0.2]), H298=(-90.466,'kcal/mol','+|-',0.1), S298=(59.464,'cal/(mol*K)','+|-',0.5), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""formic_acid""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""1 *1 C u0 p0 c0 {2,S} {3,D} {4,S}
2    O u0 p2 c0 {1,S} {5,S}
3 *2 O u0 p2 c0 {1,D}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([103.524,128.448,152.586,171.087,199.997,222.111,252.569],'J/(mol*K)'), H298=(-106.461,'kJ/mol'), S298=(371.953,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Os-CsH) + group(Os-CsH) + ring(Cyclopropane) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2     C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  *1 C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5     O u0 p2 c0 {4,S} {13,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11 *2 O u1 p2 c0 {4,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
""")])], [Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([13.614,17.846,21.434,24.32,28.647,31.786,36.685],'cal/(mol*K)','+|-',[1.3,1.6,1.6,1.6,1.5,1.4,1]), H298=(69.882,'kcal/mol','+|-',0.9), S298=(61.584,'cal/(mol*K)','+|-',0.7), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), label="""cC3H5""", comment="""Thermo library: DFT_QCI_thermo"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1    C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2    C u0 p0 c0 {1,S} {3,S} {6,S} {7,S}
3 *3 C u1 p0 c0 {1,S} {2,S} {8,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {1,S}
6    H u0 p0 c0 {2,S}
7    H u0 p0 c0 {2,S}
8    H u0 p0 c0 {3,S}
""")]), Species(label="", thermo=ThermoData(Tdata=([300,400,500,600,800,1000,1500],'K'), Cpdata=([103.524,128.448,152.586,171.087,199.997,222.111,252.569],'J/(mol*K)'), H298=(-106.461,'kJ/mol'), S298=(371.953,'J/(mol*K)'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Os-CsH) + group(Os-CsH) + ring(Cyclopropane) + radical(CCOJ)"""), molecule=[Molecule().fromAdjacencyList("""multiplicity 2
1  *3 C u0 p0 c0 {2,S} {3,S} {4,S} {6,S}
2     C u0 p0 c0 {1,S} {3,S} {7,S} {8,S}
3     C u0 p0 c0 {1,S} {2,S} {9,S} {10,S}
4  *1 C u0 p0 c0 {1,S} {5,S} {11,S} {12,S}
5     O u0 p2 c0 {4,S} {13,S}
6     H u0 p0 c0 {1,S}
7     H u0 p0 c0 {2,S}
8     H u0 p0 c0 {2,S}
9     H u0 p0 c0 {3,S}
10    H u0 p0 c0 {3,S}
11 *2 O u1 p2 c0 {4,S}
12    H u0 p0 c0 {4,S}
13    H u0 p0 c0 {5,S}
""")])]], family='R_Addition_MultipleBond', template=['CO-NdH_O', 'CsJ-CsCsH'])


    def test_cyclic_groups(self):
        groups = get_cyclic_groups(self.beta_propyl_cyclic_reaction.products[0].molecule[0])
        self.assertEqual(len(groups),1)
        self.assertIn('beta-c-3',groups.keys())
        self.assertEqual(groups['beta-c-3'],1)

        func_groups = vereecken_get_groups(self.beta_propyl_cyclic_reaction.products[0].molecule[0])
        self.assertNotIn('beta-alkyl',func_groups.keys())

    def test_improper_usage_of_group_checking(self):
        """Tests that an incorrect functional group will be found if cyclics are not checked first"""
        func_groups = vereecken_get_groups(self.beta_propyl_cyclic_reaction.products[0].molecule[0])
        #not actually supposed to be found if cyclic method called first
        self.assertIn('beta-alkyl',func_groups.keys())

    def test_formic_acid_formation_groups(self):
        groups = vereecken_get_groups(self.formic_acid_reaction.products[0].molecule[0])
        self.assertEqual(len(groups),3)
        self.assertIn('alpha-OH',groups.keys())
        self.assertIn('beta-OH',groups.keys())
        self.assertIn('beta=O',groups.keys())
        self.assertEqual(groups['beta-OH'],1)
        self.assertEqual(groups['beta=O'],1)
        self.assertEqual(groups['alpha-OH'],1)

    def test_alkyl_gamma_alcohol(self):
        groups = vereecken_get_groups(self.gamma_alcohol_beta_alkyl_rxn.products[0].molecule[0])
        self.assertEqual(len(groups),2)
        self.assertIn('alpha-alkyl-oxygenate',groups.keys())
        self.assertIn('beta-alkyl',groups.keys())
        self.assertEqual(groups['alpha-alkyl-oxygenate'],1)
    
    def test_alkyl_groups(self):
        groups = vereecken_get_groups(self.degenerate_reaction.products[0].molecule[0])
        self.assertEqual(len(groups),1)
        self.assertIn('alpha-alkyl',groups.keys())
        self.assertEqual(groups['alpha-alkyl'],1)
        
    def test_ether_groups(self):
        groups = vereecken_get_groups(self.ether_rxn.products[0].molecule[0])
        self.assertEqual(len(groups),2)
        self.assertIn('beta-C=C',groups.keys())
        self.assertIn('alpha-OR',groups.keys())
        self.assertEqual(groups['alpha-OR'],1)
        self.assertEqual(groups['beta-C=C'],1)

    def test_CO_bond_breakage_groups(self):
        with self.assertRaises(Exception):
            groups = vereecken_get_groups(self.multifunc_rxn.products[0].molecule[0])

    def test_get_arrhenius(self):
        "For ethoxyl radical decomposition"
        rate = self.sar.get_arrhenius(self.sar,self.reaction)
        self.assertAlmostEqual(rate.Ea.value,17.9,places=1) # found from atkinson2007
        self.assertAlmostEqual(rate.A.value, 1.8e13 * self.reaction.degeneracy)

    def test_formic_acid_formation_arrhenius(self):
        rate = self.sar.get_arrhenius(self.sar,self.formic_acid_reaction)
        self.assertGreater(rate.Ea.value,0)

    def test_get_arrhenius_accounts_for_degeneracy(self):
        "test CC([O])C decomposition takes into account degenerarcy"
        rate = self.sar.get_arrhenius(self.sar,self.degenerate_reaction)
        self.assertAlmostEqual(rate.A.value, 2*1.8e13)

    def test_get_rate(self):
        rate = self.sar.get_rate(self.reaction,298)
        self.assertAlmostEqual(rate, self.sar.get_arrhenius(self.sar,self.reaction).getRateCoefficient(298)) # make sure the rate is similar to in paper
        rate2 = self.sar.get_rate(self.reaction,350)
        self.assertNotAlmostEqual(rate,rate2)


if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )