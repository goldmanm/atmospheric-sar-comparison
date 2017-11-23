# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:47:20 2017

@author: mark
"""
import unittest

from atkinson2007 import sar


from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.data.kinetics.family import TemplateReaction, KineticsFamily
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.thermo import ThermoData

class SarRateEstimationTest(unittest.TestCase):

    def setUp(self):
        self.sar = sar
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


    def test_get_arrhenius(self):
        "For ethoxyl radical decomposition"
        rate = self.sar.get_arrhenius(self.sar,self.reaction)
        self.assertAlmostEqual(rate.Ea.value,18.1,places=1) # found from atkinson2007
        self.assertAlmostEqual(rate.A.value, 5e13 * self.reaction.degeneracy)

    def test_get_arrhenius_accounts_for_degeneracy(self):
        "test CC([O])C decomposition takes into account degenerarcy"
        rate = self.sar.get_arrhenius(self.sar,self.degenerate_reaction)
        self.assertAlmostEqual(rate.A.value, 2*5e13)

    def test_get_rate(self):
        rate = self.sar.get_rate(self.reaction,298)
        self.assertAlmostEqual() # make sure the rate is similar to in paper
        rate2 = self.sar.get_rate(self.reaction,350)
        self.assertNotAlmostEqual(rate,rate2)

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )