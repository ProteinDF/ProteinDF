#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse

from bratomgroup import *


class BrPdb(object):
    def __init__(self):
        """
        create empty PDB object
        """
        self.__data = []

    def get_number_of_items(self):
        return len(self.__data)


    def get_name(self, index):
        answer = None
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('name', None)
        return answer


    def set_name(self, index, name):
        if ((0 <= index) and (index < self.get_number_of_items())):
            self.__data[index]['name'] = name


    def get_element(self, index):
        answer = None
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('element', None)
        return answer


    def set_element(self, index, element):
        if ((0 <= index) and (index < self.get_number_of_items())):
            self.__data[index]['element'] = element


    def get_posision(self, index):
        answer = None
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('coord', None)
        return answer


    def set_position(self, index, position):
        if ((0 <= index) and (index < self.get_number_of_items())):
            self.__data[index]['coord'] = position


    def get_charge(self, index):
        answer = 0
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('charge', 0)
        return answer


    def set_charge(self, index, charge):
        if ((0 <= index) and (index < self.get_number_of_items())):
            self.__data[index]['charge'] = charge


    def get_occupancy(self, index):
        answer = None
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('occupancy', None)
        return answer


    def get_temp_factor(self, index):
        answer = None
        if ((0 <= index) and (index < self.get_number_of_items())):
            answer = self.__data[index].get('temp_factor', None)
        return answer
        

    def set_temp_factor(self, serial, temp_factor):
        serial = int(serial)
        self.__data[serial]['temp_factor'] = temp_factor


    def renumber(self):
        for index in range(len(self.__data)):
            self.__data[index]['serial'] = index +1

    def sort_by_serial(self):
        self.__data.sort(cmp = lambda x, y: cmp(int(x['serial']), int(y['serial'])))

    def load(self, file_path):
        if (os.path.isfile(file_path) != True):
            return

        fin = open(file_path, "r")
        while True:
            line = fin.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')
            
            #if (len(line) != 80):
            #    continue

            record_name = line[0:6]
            if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):
                if (len(line) < 80):
                    line = line + (' ' * (80 - len(line)))

                serial = int(line[6:11])
                name = line[12:16]
                alt_loc = line[16]
                res_name = line[17:20]
                chain_id = line[21]
                res_seq = line[22:26]
                i_code = line[26]
                coord_x = line[30:38]
                coord_y = line[38:46]
                coord_z = line[46:54]
                occupancy = line[54:60].strip()
                temp_factor = line[60:66].strip()
                element = line[76:78].strip()
                charge = line[78:80].strip()
            
                item = {}
                item['serial'] = serial
                item['record_name'] = record_name
                item['name'] = name
                item['alt_loc'] = alt_loc
                item['res_name'] = res_name
                item['chain_id'] = chain_id
                item['res_seq'] = int(res_seq)
                item['i_code'] = i_code
                item['coord'] = [float(coord_x), float(coord_y), float(coord_z)]

                if (len(occupancy) != 0):
                    item['occupancy'] = float(occupancy)
                else:
                    item['occupancy'] = 1.0
                if (len(temp_factor) != 0):
                    item['temp_factor'] = float(temp_factor)
                else:
                    item['temp_factor'] = 0.0

                if (len(element) != 0):
                    item['element'] = element
                else:
                    if (name == ' MG '):
                        element = 'Mg'
                    elif (name == 'FE  '):
                        element = 'Fe'
                    else:
                        element = name[1]
                    item['element'] = element

                if (len(charge) != 0):
                    item['charge'] = charge
                else:
                    item['charge'] = "  "

                self.__data.append(item)
                continue
            elif (record_name == 'TER   '):
                serial = line[6:11]
                resname = line[17:20]
                chain_id = line[21]
                res_seq = line[22:26]
                i_code = line[26]
                item = {}
                item['serial'] = int(serial)
                item['record_name'] = record_name
                item['res_name'] = res_name
                item['chain_id'] = chain_id
                item['res_seq'] = res_seq
                item['i_code'] = i_code

                self.__data.append(item)
                continue


    def get_atom_group(self):
        root = BrAtomGroup()
        for index in range(len(self.__data)):
            item = self.__data[index]
            record_name = item['record_name']
            serial = item['serial']
            if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):
                name = item['name'].strip()
                alt_loc = item['alt_loc']
                res_name = item['res_name']
                chain_id = item['chain_id']
                res_seq = str(item['res_seq'])
                i_code = item['i_code']
                coord = item['coord']
                occupancy = item.get('occupancy', 1.0)
                temp_factor = item.get('temp_factor', 0.0)
                element = item.get('element', 'X')
                charge = item.get('charge', 0.0)

                if (chain_id == " "):
                    chain_id = "_"

                if (root.has_group(chain_id) == False):
                    chain = BrAtomGroup()
                    chain.set_name(chain_id)
                    root.set_group(chain_id, chain)
                if (root[chain_id].has_group(res_seq) == False):
                    residue = BrAtomGroup()
                    residue.set_name(res_name)
                    root[chain_id].set_group(res_seq, residue)
                atom = BrAtom()
                atom.set_element(element)
                atom.set_position(BrPosition(coord))
                atom.set_name(name)
                root[chain_id][res_seq].set_atom(serial, atom)
        return root

    def set_by_atomgroup(self, atomgroup):
        assert(isinstance(atomgroup, BrAtomGroup) == True)

        self.__data = []
        item = {}
        item["record_name"] = "ATOM  "
        item["alt_loc"] = " "
        item["i_code"] = " "
        for chain_id in atomgroup.get_group_list():
            if (chain_id != "_"):
                item['chain_id'] = chain_id
            else:
                item['chain_id'] = " "
            chain = atomgroup.get_group(chain_id)

            for res_seq in chain.get_group_list():
                residue = chain.get_group(res_seq)
                item["res_seq"] = res_seq
                item["res_name"] = residue.get_name()

                for serial in residue.get_atom_list():
                    atom = residue.get_atom(serial)
                    if (atom != None):
                        item["serial"] = serial
                        name = atom.get_name()
                        if (len(name) != 4):
                            name = self.__match_name_table(atom.get_name(), atom.get_symbol())
                        item["name"] = name
                        item["coord"] = atom.get_position().get_list()
                        self.__data.append(copy.deepcopy(item))
        self.sort_by_serial()
        
    def __str__(self):
        occupancy = 1.0
        output = ""
        #for serial, item in self.__data.items():
        for index in range(len(self.__data)):
            item = self.__data[index]
            record_name = item['record_name']
            serial = int(item['serial'])
            if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):
                name = item['name']
                alt_loc = item['alt_loc']
                res_name = item['res_name']
                chain_id = item['chain_id']
                res_seq = int(item['res_seq'])
                i_code = item['i_code']
                coord = item['coord']
                occupancy = item.setdefault('occupancy', 1.0)
                temp_factor = item.setdefault('temp_factor', 1.0)
                element = item.setdefault('element', '  ')
                charge = item.setdefault('charge', '  ')
                line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % (
                    serial, name, alt_loc,
                    res_name, chain_id, res_seq, i_code,
                    coord[0], coord[1], coord[2],
                    occupancy, temp_factor, element.upper(), charge)
                output += line
            elif (record_name == 'TER   '):
                line = "TER   %5d      %3s %c%4d%c\n" % (serial, res_name, chain_id, res_seq, i_code)
                output += line

        return output

    def __match_name_table(self, name, symbol):
        assert(isinstance(name, str) == True)
        assert(isinstance(symbol, str) == True)
        capital_name = name.upper()
        capital_symbol = symbol.upper()
        if (capital_symbol == "NA"):
            name = "NA  "
        elif (capital_symbol == "CL"):
            name = "CL  "
        elif (capital_name[0] == capital_symbol[0]):
            name = " %s" % (name)
        name += " " * (4 - len(name))
        return name


    def __iter__(self):
        for item in self.__data:
            yield item
        

def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] PDB_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output_path",
                      help = "PDB output file", metavar = "FILE")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    file_path = args[0]
    verbose = opts.verbose

    # 
    pdb_obj = BrPdb(file_path)
    print(pdb_obj)

    # end

if __name__ == '__main__':
    main()

