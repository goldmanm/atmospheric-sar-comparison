# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 15:39:29 2017

@author: mark
"""

import unittest
from sar import SAR
from rmgpy.molecule import Molecule

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