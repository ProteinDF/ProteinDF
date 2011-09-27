#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import copy
import optparse
import profile

from bratomgroup import *


class BrComponentsCif(object):
    def __init__(self):
        """
        """
        self.__data = {}

    def load(self, file_path):
        if (os.path.isfile(file_path) != True):
            return

        # for data_set
        data = ""
        loop_head = []
        loop_mode = False
        stored_line = ""

        fin = open(file_path, "r")
        while True:
            line = fin.readline()
            if (len(line) == 0):
                break
            line = line.rstrip("\n")
            #print(line)
            
            if (line[0] == "#"):
                if (loop_mode == True):
                    loop_mode = False
                    stored_line = ""
                    loop_head = []
                    re_loop_body = None
                continue

            # data_* section
            if (line[0:5] == "data_"):
                data = line[5:]
                loop_mode = False
                re_loop_body = None
                continue

            if (loop_mode == False):
                if (line[0] == "_"):
                    word_list = self.get_word_list(line)
                    while (len(word_list) < 2):
                        next_line = fin.readline()
                        next_line = next_line.rstrip("\n")
                        #print(next_line)
                        # special operation for ";"
                        if (next_line[0] == ';'):
                            while (True):
                                cont_line = fin.readline()
                                cont_line = cont_line.rstrip()
                                #print(cont_line)
                                next_line += cont_line
                                if (cont_line == ";"):
                                    next_line = next_line[1:-1]
                                    break
                        line = line + next_line
                        word_list = self.get_word_list(line)
                    #print(line, word_list)
                    keyword = word_list[0]
                    keyword = keyword[1:] # remove "_"
                    key_list = keyword.split(".", 1)
                    key1 = key_list[0]
                    key2 = key_list[1]
                    value = word_list[1]
                    self.__data.setdefault(data, {})
                    self.__data[data].setdefault(key1, {})
                    self.__data[data][key1][key2] = value
                elif (line == "loop_"):
                    loop_mode = True
                    loop_head = []
                else:
                    # something wrong
                    print("unsupport format: %s" % (line))

            else:
                #loop mode
                if (len(stored_line) > 0):
                    line = stored_line + line
                    stored_line = ""

                if (line[0] == "_"):
                    # loop header
                    line = line.rstrip()
                    loop_head.append(line[1:])
                    continue
                else:
                    # loop body
                    line = stored_line + line
                    word_list = self.get_word_list(line)
                    if (len(word_list) < len(loop_head)):
                        stored_line = line
                        continue
                    else:
                        item = {}
                        key1 = ""
                        for index, head_str in enumerate(loop_head):
                            key_list = head_str.split(".", 1)
                            key1 = key_list[0]
                            key2 = key_list[1]
                            item[key2] = word_list[index]
                        self.__data.setdefault(data, {})
                        self.__data[data].setdefault(key1, [])
                        self.__data[data][key1].append(copy.deepcopy(item))
                    

    def get_word_list(self, line):
        """
        スペースで区切られた文字列(単語)のリストを返す。
        ダブルクォートで囲まれている単語はきちんとつなげる。
        """
        answer = []
        parts = line.split()

        index = 0
        max_index = len(parts)
        continue_mode = False
        stored_line = ""
        bracket_char = ""
        while (index < max_index):
            value = parts[index]
            index += 1
            
            if (continue_mode == False):
                check_char = value[0]
                if (check_char != '"'):
                    answer.append(value)
                    continue
                else:
                    continue_mode = True
                    bracket_char = value[0]
                    # go below

            # continue_mode
            stored_line += " " + value
            stored_line = stored_line.strip()
            if (stored_line[-1] == bracket_char):
                continue_mode = False
                answer.append(stored_line[1:-1])
                stored_line = ""

        return answer


    def __str__(self):
        output = ""
        for data_name in self.__data.keys():
            output += "data: %s\n" % (data_name)
            for category, item in self.__data[data_name].items():
                if (isinstance(item, list) == True):
                    for itr, subcategory in enumerate(item):
                        for name, value in subcategory.items():
                            output += " [%s][%d][%s] = [%s]\n" % (category, itr, name, value)
                elif (isinstance(item, dict) == True):
                    for name, value in subcategory.items():
                        output += " [%s][%s] = [%s]\n" % (category, name, value)
                else:
                    output += " [%s] = [%s]\n" % (category, item)
            output += "\n"

        return output



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
    cif_obj = BrComponentsCif()
    cif_obj.load(file_path)
    print(cif_obj)
    
    # end

if __name__ == '__main__':
    main()
    #profile.run("main()")

