# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 09:14:26 2017

@author: mark
"""

from sar import AlkoxyDecompositionSAR
from sarDatabase import SARDatabase

sar = AlkoxyDecompositionSAR(label='Atkinson2007',
                         valid_atom_types = ['Cs','CO',''],
                         temperature_range = (244,312.5),
                         pressure_range = (1e5,1e5),
                )

sar.a_value_tree = SARDatabase().load('/home/mark/workspace/alkoxy-sar/sar_data/atkinson2007.py')

def atkinson_get_arrhenius(self,reaction):
    """
    Estimates the reaction rate for alkoxy decomposition given
    an RMG reaction from R_Addition_multiplebond family
    using species with proper labeling.
    
    The reaction object will be in the oposite direction with the product being
    alkoxy radical

    Returns an RMG kineticsdata object

    The method originates from Atkinson 2007
    """
    from rmgpy.kinetics.arrhenius import Arrhenius

    A = 5e13 #s-1
    b = 0.40
    #get radical species for lookup of a value
    radical_decomp_product = None
    for s in reaction.reactants:
        mol = s.molecule[0]
        labeled_atoms = mol.getLabeledAtoms()
        try:
            radical_atom = labeled_atoms['*3']
            radical_decomp_product = mol
            break
        except KeyError:
            pass
    if radical_decomp_product is None:
        if len(labeled_atoms) == 0:
            raise TypeError("atoms not labeled for this reaction")
        raise TypeError('Reaction object does not have specified product. RXN: {}'.format(repr(reaction)))
    tree_entry = self.a_value_tree.descendTree(radical_decomp_product,labeled_atoms)
    a = tree_entry.data
    if a is None:
        raise TypeError('The node {} had no data'.format(tree_entry))

    #get heat of reaction
    deltaH = reaction.getEnthalpyOfReaction(298) # j/mol
    deltaH /=4184 #kcal/mol

    Ea = a - b*deltaH

    # get deneracy
    A *= reaction.degeneracy

    return Arrhenius(A=(A,'s^-1'),Ea=(Ea,'kcal/mol'))

sar.get_arrhenius = atkinson_get_arrhenius