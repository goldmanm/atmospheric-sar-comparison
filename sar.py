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
        self._check_proper_attribute_types(self)

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

    def get_rate(self, molecules,temperature, pressure, **kwargs):
        """
        input list of molecular graph object from RMG-mol with the reacting atoms labeled, 
        temperature (kevlin), pressure (pascals), and other auxilliary information

        outputs the rate constant (float) at those conditions
        """
        raise NotImplementedError('Please overwrite get_rate method')

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
        alkohols and ketones but not both together in the same molecule. 
        """
        if len(self.valid_atom_types) == 0:
            raise ValueError('Cannot test non_valid_types of molecule {0} if self.valid_atom_types is empty'.format(molecule))
        non_valid_types = set()
        for atom in molecule.atoms:
            if atom.type not in self.valid_atom_types:
                non_valid_types.add(atom.type)
        return non_valid_types

class AlkoxySAR(SAR):
    """
    Inherited class from SAR for specifically Alkoxy radical isomerizations.
    """