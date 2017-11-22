# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 15:39:29 2017

@author: mark
"""

import unittest
from sar import SAR
from rmgpy.molecule import Molecule
from sarDatabase import SARDatabase
from sar import get_labeled_atoms_alkoxy_decomposition
from sarInstances import sar_atkinson

from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.data.kinetics.family import TemplateReaction, KineticsFamily
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.thermo import ThermoData

class SarTest(unittest.TestCase):
    
    def setUp(self):
        molecule_adj_lists = ["""
        
        """,
        """
        
        """]
        
        #self.molecules = [Molecule.fromAdjacencyList(x) for x in molecule_adj_list]


    def test_initialization(self):
        self.sar = SAR(label = 'testaSARous',
                       valid_atom_types = ['Cs','CO'],)
        self.assertEqual(self.sar.label, 'testaSARous')
        self.assertIn('CO', self.valid_atom_types)


if __name__ =='__main__':
    unittest.main()
    
    
class SarDatabaseTest(unittest.TestCase):

    def setUp(self):
        self.database = SARDatabase().load('sar_data/Atkinson2007.py')
        self.methyl = Molecule().fromAdjacencyList("""
multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
""")
        self.nbutyl = Molecule().fromAdjacencyList("""
multiplicity 2
1 *3 C u1 p0 c0 {2,S} {3,S} {4,S}
2  H u0 p0 c0 {1,S}
3  H u0 p0 c0 {1,S}
4  C u0 p0 c0 {1,S} {5,S} {7,S} {8,S}
5  C u0 p0 c0 {4,S} {6,S} {9,S} {10,S}
6  C u0 p0 c0 {5,S} {11,S} {12,S} {13,S}
7  H u0 p0 c0 {4,S}
8  H u0 p0 c0 {4,S}
9  H u0 p0 c0 {5,S}
10 H u0 p0 c0 {5,S}
11 H u0 p0 c0 {6,S}
12 H u0 p0 c0 {6,S}
13 H u0 p0 c0 {6,S}
""")

        self.nitrogen = Molecule().fromAdjacencyList("""
multiplicity 2
1 C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}
2 *3 N u1 p1 c0 {1,S} {3,S}
3 C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {1,S}
7 H u0 p0 c0 {3,S}
8 H u0 p0 c0 {3,S}
9 H u0 p0 c0 {3,S}

""")


        self.ethane = Molecule().fromAdjacencyList("""
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 *3 C u0 p0 c0 {1,S} {6,S} {7,S} {8,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
7 H u0 p0 c0 {2,S}
8 H u0 p0 c0 {2,S}
""")


    def test_lookup_molecule_entry(self):
        atoms = get_labeled_atoms_alkoxy_decomposition(self.methyl)
        self.assertNotEqual(len(atoms),0)
        entry = self.database.descendTree(self.methyl,atoms)
        self.assertEqual(entry.label,'methyl')
        self.assertEqual(entry.data,12.9)

    def test_lookup_group_entry(self):
        atoms = get_labeled_atoms_alkoxy_decomposition(self.nbutyl)
        self.assertNotEqual(len(atoms),0)
        entry = self.database.descendTree(self.nbutyl,atoms)
        self.assertEqual(entry.label,'primary_C')
        self.assertEqual(entry.data,9.5)

    def test_lookup_head_entry_no_data(self):
        atoms = get_labeled_atoms_alkoxy_decomposition(self.nitrogen)
        self.assertNotEqual(len(atoms),0)
        entry = self.database.descendTree(self.nitrogen,atoms)
        self.assertEqual(entry.label,'parent')
        self.assertIsNone(entry.data)

    def test_lookup_incorrect_entry(self):
        atoms = get_labeled_atoms_alkoxy_decomposition(self.ethane)
        self.assertNotEqual(len(atoms),0)
        entry = self.database.descendTree(self.ethane,atoms)
        self.assertIsNone(entry)


class SarRateEstimationTest(unittest.TestCase):

    def setUp(self):
        self.sar = sar_atkinson
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
    def test_lookup_molecule_entry(self):
        rate = self.sar.get_arrhenius(self.sar,self.reaction)
        self.assertAlmostEqual(rate.Ea.value,18.1,places=1) # found from atkinson2007
        self.assertFail()

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
