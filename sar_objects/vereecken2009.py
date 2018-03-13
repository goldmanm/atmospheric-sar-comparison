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
    #print alkoxy
    #print groups
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

def get_cyclic_groups(alkoxy_mol):
    """
    Given a molecule object, returns all the groups vereecken used
    to represent the molecules as a dictionary with keys being strings
    and values being the frequency of groups being formed.

    This does not account for breaking of rings during decomposition
    """
    labeled_atoms = alkoxy_mol.getLabeledAtoms()
    alpha_cycles = alkoxy_mol.getAllCycles(labeled_atoms['*1'])
    beta_cycles = alkoxy_mol.getAllCycles(labeled_atoms['*3']) # returns duplicate of each cycle
    cycle_groups = {}
    groups_used = []
    for cycle in alpha_cycles + beta_cycles:
        length = len(cycle)
        alpha = False
        beta = False
        group_to_add = ""
        if any([a.label == '*1' for a in cycle]):
            alpha=True
        if any([a.label == '*3' for a in cycle]):
            beta=True
        if alpha and beta:
            print('cycle for molecule {} contains a ring in both alpha and beta'.format(alkoxy_mol))
        elif alpha and length >3 and length <7:
            group_to_add = 'alpha-c-{}'.format(length)
        elif beta and length > 2 and length < 7:
            group_to_add = 'beta-c-{}'.format(length)

        if group_to_add != '':
            try:
                # 0.5 returned since each cycle is duplicated
                cycle_groups[group_to_add] += 0.5
            except KeyError:
                cycle_groups[group_to_add] = 0.5
            groups_used.append(cycle)

    # now remove the bonds between the counted cycles to prevent counting
    # them as functional groups as well
    for cycle in groups_used:
        for atom_index in range(len(cycle)):
            if alkoxy_mol.hasEdge(cycle[atom_index],cycle[atom_index-1]):
                edge = alkoxy_mol.getEdge(cycle[atom_index],cycle[atom_index-1])
                alkoxy_mol.removeEdge(edge)
    return cycle_groups

def vereecken_get_groups(alkoxy_mol):
    # make a copy for analysis, since get_cyclic_groups can break bonds
    # to prevent double counting with functional groups
    alkoxy = alkoxy_mol.copy(deep=True)
    cyclic_groups = get_cyclic_groups(alkoxy)
    groups = get_functional_groups(alkoxy)
    groups.update(cyclic_groups)
    return groups

def get_functional_groups(alkoxy_mol):
    """
    given a molecule object `alkoxy_mol`. This method returns
    a dictionary of groups used in the Vereecken SAR with the
    key being the group and the value being the number of occurances
    it has.
    """
    #print 'getting groups from {}'.format(alkoxy_mol.toSMILES())
    alkoxy_mol.assignAtomIDs()
    labeled_atoms = alkoxy_mol.getLabeledAtoms()
    assert labeled_atoms['*1'].symbol == 'C'
    assert labeled_atoms['*3'].symbol == 'C', alkoxy_mol.toAdjacencyList() + str(labeled_atoms)
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
                    found_valid_group = False
                    for nnn,bond2 in nn.bonds.items():
                        if nnn.label == '':
                            if nnn.isOxygen():
                                if bond2.isDouble():
                                    found_valid_group=True
                                else:
                                    for nnnn in nnn.bonds.keys():
                                        if nnnn.id != nn.id:
                                            if nnnn.isOxygen() or nnnn.isCarbon():
                                                found_valid_group = True
                    if found_valid_group:
                        addon = '-alkyl-oxygenate'
                    else:
                        addon = '-alkyl'
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
            elif nn.symbol != 'H':
                print('non-valid atom type found for atom {} on atom {} with bonds {}'.format(nn,atom, atom.bonds))
    if '-alkyl-oxygenate' in functional_groups.keys() and len(functional_groups)>1:
        if '-alkyl' in functional_groups.keys():
            functional_groups['-alkyl'] += functional_groups['-alkyl-oxygenate']
        else:
            functional_groups['-alkyl'] = functional_groups['-alkyl-oxygenate']
        del functional_groups['-alkyl-oxygenate']
    return functional_groups

sar.get_arrhenius = vereecken_get_arrhenius