import csv

from paramagpy import dataparse


class PCSReader:
    # __result = []
    __values = dataparse.DataContainer(dtype='PCS')

    def __init__(self, file_path, predicate=None):
        self.__file = open(file_path)
        reader = csv.DictReader(self.__file, dialect='excel-tab')
        for row in reader:
            # self.__result.append(PCSAtom(row['Residue'], row['Residue_Nr'], row['Group'], row['Atom'],
            # row['PCS_exp']))
            if predicate is None or predicate(row):
                self.__values[int(row['Residue_Nr']), row['Atom']] = float(row['PCS_exp']), float(0)

    def get_result(self):
        # return self.__result
        return self.__values
