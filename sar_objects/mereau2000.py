# -*- coding: utf-8 -*-

from sar import AlkoxyDecompositionSAR
from sarDatabase import SARDatabase
import pandas as pd
from orlando2003 import get_IP_of_reaction, get_ionization_potential

sar = AlkoxyDecompositionSAR(label=u'Méreau2003',
                         valid_atom_types = ['Cs','CO',''],
                         temperature_range = (244,312.5),
                         pressure_range = (1e5,1e5),
                )

def mereau_get_arrhenius(self,reaction):
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

    A = 1e14 #s-1
    IP = get_IP_of_reaction(reaction)

    n_H = get_n_H(reaction.products[0].molecule[0])
    Ea = 2.5 * IP +2.1*n_H - 10.4

    # get deneracy
    A *= reaction.degeneracy

    return Arrhenius(A=(A,'s^-1'),Ea=(Ea,'kcal/mol'))


sar.get_arrhenius = mereau_get_arrhenius

def get_n_H(alkoxy_rad):
    """
    outputs the number of hydrogens attached to the carbon
    closest to the alkoxy radical group.
    """
    labeled_atoms = alkoxy_rad.getLabeledAtoms()
    alpha_carbon_to_oxy_rad = labeled_atoms['*1']
    n_hydrogens = 0
    for atom in alpha_carbon_to_oxy_rad.bonds.keys():
        if atom.isHydrogen():
            n_hydrogens += 1
    return n_hydrogens