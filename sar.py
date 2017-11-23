# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 15:39:31 2017

@author: mark
"""

class SAR:
    """
    attributes:
    
    label               string identifying the sar
    description         string with detailed information about the SAR (source, etc)
    valid_atom_types    list of strings describing the valid atom types of the sar
    temperature_range   length 2 tuple of doubles with the minimum and maximum valid temperatures (K)
    pressure_range      length 2 tuple of doubles with min and max pressures for model (Pa)
    
    methods:
    
    get_rate            
    is_valid_reactant   
    non_valid_types     
    """
    def __init__(self,label = '',
                 valid_atom_types = None,
                 temperature_range=None,
                 pressure_range = None,
                 description = ''
                 ):
        self.label = label
        self.valid_atom_types = valid_atom_types or []
        self.temperature_range = temperature_range or (200,1000)
        self.pressure_range = pressure_range or (1e4,1e6)
        self.description = description
        self._check_proper_attribute_types()

    def _check_proper_attribute_types(self):
        """
        checks that attributes are in the propper type. Meant to catch most, but
        not all errors
        """
        assert isinstance(self.label,str)
        assert isinstance(self.valid_atom_types, list)
        if self.valid_atom_types:
            assert isinstance(self.valid_atom_types[0], str)
        assert len(self.temperature_range) == 2
        assert len(self.pressure_range) == 2

    def get_rate(self, reaction,temperature, pressure=1e5, **kwargs):
        """
        input list of molecular graph object from RMG-mol with the reacting atoms labeled, 
        temperature (kevlin), pressure (pascals), and other auxilliary information

        outputs the rate constant (float) at those conditions
        """
        rate_expression = self.get_arrhenius(self, reaction, **kwargs)
        return rate_expression.getRateCoefficient(temperature, pressure)

    def get_arrhenius(self, reaction, **kwargs):
        """
        input list of molecular graph objects from RMG-mol that represent the reactant
        and any other auxilliary arguments needed for this method.

        outputs RMG-kinetics Arrhenius object for the given reaction.
        """
        raise NotImplementedError('Please overwrite get_rate or get_arrhenius method for SAR {}'.format(self.label))

    def is_valid_reactant(self,molecule):
        """
        checks that all the atom_types in the molecule are in the valid_atom_types
        list

        returns a tuple with a boolean (true if the reactant is valid) followed by 
        a set of all the non-valid types()
        """
        if self.non_valid_types(self,molecule):
            return False
        return True

    def non_valid_types(self,molecule):
        """
        checks that all the atom_types in the molecule are in the valid_atom_types
        list

        returns a set of all atom types in molecule that are not regressed in the SAR

        note: this may need to be further refined since a model may be able to support
        alcohols and ketones but not both together in the same molecule. also does not
        check cyclics.
        """
        if len(self.valid_atom_types) == 0:
            raise ValueError('Cannot test non_valid_types of molecule {0} if self.valid_atom_types is empty'.format(molecule))
        non_valid_types = set()
        for atom in molecule.atoms:
            if atom.type not in self.valid_atom_types:
                non_valid_types.add(atom.type)
        return non_valid_types

class AlkoxyDecompositionSAR(SAR):
    """
    Inherited class from SAR for specifically Alkoxy radical decompositions
    """

    # helper methods for getting rates from various instances

def get_labeled_atoms_alkoxy_decomposition(molecule):
    """
    input is a reaction of molecule objects from RMG-mol with proper labeling

    output is a dictionary with keys 'radical', 'carbonyl oxygen' and 'carbonyl carbon'
    """
    output = {}
    for atom in molecule.atoms:
        if atom.label == '*3':
            output['radical'] = atom
        elif atom.label == '*1':
            output['carbonyl carbon'] = atom
        elif atom.label == '*2':
            output['carbonyl oxygen'] = atom
    return output

class AlkoxyIsomerizationSAR(SAR):
    """
    Inherited class from SAR for specifically Alkoxy radical decompositions
    
    Not complete
    """
    def __init__(self,label = '',
                 valid_atom_types = None,
                 temperature_range=None,
                 pressure_range = None,
                 description = '',
                 valid_distances = None
                 ):
        self._SAR__init__(label, valid_atom_types, temperature_range, pressure_range,description)
        self.valid_distances = valid_distances or []

    # helper methods for getting rates from various instances
    def get_labeled_atoms(molecule):
        """
        input is a molecular graph from RMG-mol with proper labeling

        output is a dictionary with keys 'radical', 'hydrogen' and 'alkane'
        """
        output = {}
        for atom in molecule.atoms:
            if atom.label == '*1':
                output['radical'] = atom
            elif atom.label == '*2':
                output['alkane'] = atom
            elif atom.label == '*3':
                output['hydrogen'] = atom
        return output

    def get_ts_ring_size(molecule):
        """
        input is a molecular graph from RMG-mol with labeled atoms

        output is an integer describing the number of atoms in final transitionstate
        """
        pass