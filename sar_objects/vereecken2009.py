# -*- coding: utf-8 -*-


from sar import AlkoxyDecompositionSAR
from sarDatabase import SARDatabase

sar = AlkoxyDecompositionSAR(label='Vereeken2007',
                         valid_atom_types = ['Cs','CO',''],
                         temperature_range = (0,1000),
                         pressure_range = (1e5,1e5),
                )

def vereecken_get_arrhenius(self,reaction):
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

    A = 1.8e13 #s-1
    n = 1.7
    T0 = 298

    #get alkoxy
    alkoxy = reaction.products[0].molecule[0]
    groups = vereecken_get_groups(alkoxy)
    Eb = calculate_Eb(groups)

    # get deneracy
    A *= reaction.degeneracy

    return Arrhenius(A=(A,'s^-1'),n= n, T0=(T0,'K'),Ea=(Eb,'kcal/mol'))

import pandas as pd
lookup_table = pd.read_csv('/home/mark/workspace/alkoxy-sar/sar_data/vereecken2009.csv',index_col='group')['effect']

def calculate_Eb(groups):
    base = 17.9
    for label, number in groups.items():
        base += lookup_table[label] * number
    if base >= 7:
        return base
    # we need to add modification
    import math
    new_base = 19*math.exp(-(base - 22)**2/225)
    return new_base

def vereecken_get_groups(alkoxy_mol):
    """
    given a molecule object `alkoxy_mol`. This method returns
    a dictionary of groups used in the Vereecken SAR with the
    key being the group and the value being the number of occurances
    it has.
    """
    alkoxy_mol.assignAtomIDs()
    labeled_atoms = alkoxy_mol.getLabeledAtoms()
    assert labeled_atoms['*1'].symbol == 'C'
    assert labeled_atoms['*3'].symbol == 'C'
    alpha_groups = get_atom_groups(labeled_atoms['*1'])
    beta_groups = get_atom_groups(labeled_atoms['*3'])
    
    # find cyclic groups here (after project finished)
    
    all_groups = {}
    for label, num in alpha_groups.items():
        all_groups['alpha{}'.format(label)] = num
    for label, num in beta_groups.items():
        all_groups['beta{}'.format(label)] = num
    return all_groups

def get_atom_groups(atom):
    """
    returns a dictionary of functional groups attached to the
    specified atom `atom` for molecule `molecule`, which contains
    atom ID numbers. Ignores neighbors with labels.
    
    Does not get cyclic values.
    """
    functional_groups = {}
    for nn, bond in atom.bonds.items():
        if nn.label == '':
            addon = ""
            # only look at unlabeled atoms
            if nn.isCarbon():
                if bond.isDouble():
                    addon = '=C'
                elif nn.atomType.label == 'Cd':
                    addon = '-C=C'
                elif any([nnn.isOxygen() for nnn in nn.bonds.keys() if nnn.label=='']):
                    addon = '-alkyl-oxygenate'
                else:
                    addon = '-alkyl'

            elif nn.isOxygen():
                if len(nn.bonds) ==2:
                    nnn = [nnn for nnn in nn.bonds.keys() if nnn.id != atom.id][0]
                    if nnn.isHydrogen():
                        addon = '-OH'
                    elif nnn.isCarbon():
                        addon = '-OR'
                    elif nnn.isNitrogen():
                        bond_num = len(nnn.bonds)
                        if bond_num == 3:
                            addon = '-ONO2'
                        elif bond_num == 2:
                            addon = '-ONO'
                    elif nnn.isOxygen():
                        nnnn = [nnnn for nnnn in nnn.bonds.keys() if nnnn.id != nn.id][0]
                        if nnnn.isHydrogen():
                            addon='-OOH'
                        elif nnnn.isCarbon():
                            addon = '-OOR'
                elif bond.isDouble():
                    addon = '=O'

            elif nn.isNitrogen():
                bond_num = len(nn.bonds)
                if bond_num == 3:
                    addon='-NO2'
                elif bond_num == 2:
                    addon = '-NO'
            if addon != '':
                try:
                    functional_groups[addon] += 1
                except KeyError:
                    functional_groups[addon] = 1
            else:
                print('non-valid atom type found on atom {} with bonds {}'.format(atom, atom.bonds))
    return functional_groups

sar.get_arrhenius = vereecken_get_arrhenius