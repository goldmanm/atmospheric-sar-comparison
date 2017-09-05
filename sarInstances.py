# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 09:14:26 2017

@author: mark
"""

from sar import AlkoxySAR

sar_instance = AlkoxySAR(label='Atkinson2007',
                         valid_atom_types = ['Cs','CO',''],
                         temperature_range = (244,312.5),
                         pressure_range = (1e5,1e5),
                         valid_distances = [6]
                )

def atkinson_get_arrhenius(self, molecules):
    molecule = molecules[0]
    