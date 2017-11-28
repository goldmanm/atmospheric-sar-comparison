# -*- coding: utf-8 -*-


from sar import AlkoxyDecompositionSAR
from sarDatabase import SARDatabase
import pandas as pd
sar = AlkoxyDecompositionSAR(label='Orlando2003',
                         valid_atom_types = ['Cs','CO',''],
                         temperature_range = (244,312.5),
                         pressure_range = (1e5,1e5),
                )

def orlando_get_arrhenius(self,reaction):
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

    #get heat of reaction
    deltaH = reaction.getEnthalpyOfReaction(298) # j/mol
    deltaH /=4184 #kcal/mol

    Ea = 2.4 * IP - 0.58*deltaH -11.8

    # get deneracy
    A *= reaction.degeneracy

    return Arrhenius(A=(A,'s^-1'),Ea=(Ea,'kcal/mol'))

def get_IP_of_reaction(reaction):
    for spec in reaction.reactants:
        labeled_atoms = spec.molecule[0].getLabeledAtoms()
        try:
            labeled_atoms['*3']
            radical = spec
        except KeyError:
            pass
    # now check for ionization potential matches for each structure
    for mol in radical.molecule:
        try:
            return get_ionization_potential(mol)
        except KeyError:
            pass
    raise Exception("No reaction found for species {} in reaction {}".format(radical,reaction))

# create a molecule to IP dictionary for use in get_ionization_potentials
from rmgpy.molecule import Molecule
nist_data_table = pd.read_csv('/home/mark/workspace/alkoxy-sar/sar_data/IP.csv')
nist_ionization_potentials = {}
for i in nist_data_table.index:
    nist_ionization_potentials[Molecule(SMILES=nist_data_table.loc[i,'smiles'])] =\
                                nist_data_table.loc[i,'IP']

def get_ionization_potential(radical):
    """
    given a molecule object, finds IP value from the table lookup
    NIST values.

    Return a float in units of eV
    """
    try:
        return nist_ionization_potentials[radical]
    except KeyError:
        raise KeyError('Radical {} not found in list of ionization potentials.'.format(radical.toSMILES()))

sar.get_arrhenius = orlando_get_arrhenius