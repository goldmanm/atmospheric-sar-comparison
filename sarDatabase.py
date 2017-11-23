# -*- coding: utf-8 -*-

from rmgpy.data.base import Database, Entry
from rmgpy.molecule import Molecule
from rmgpy.molecule.group import Group

class SARDatabase(Database):
    """
    This class has a tree structure that extends from Database, which is used
    for finding the proper node.
    """

    def loadEntry(self, label, molecule=None, group=None, data=None):
        """
        Load an entry from the forbidden structures database. This method is
        automatically called during loading of the forbidden structures 
        database.
        """
        assert molecule is not None or group is not None
        assert not (molecule is not None and group is not None)
        if molecule is not None:
            item = Molecule().fromAdjacencyList(molecule)
        elif group is not None:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            label = label,
            item = item,
            data = data,
        )