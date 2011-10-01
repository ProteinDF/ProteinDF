#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import time
import copy

from dfdata import *

class DfOutput(object):
    """
    parse and analize ProteinDF output file(eg. fl_Out_Std).
    """
    #
    # compiled regular expressions =====================================
    #
    # step
    __re_step = re.compile("^\s*>>>>\s*(\S.*\S)\s*\(\S.+\)$")
    
    # input section
    __re_input_data_start_section = re.compile(">>>>\s*INPUT\s*DATA")
    __re_input_data_end_section = re.compile("^\s*>>>>\s*(INTEGRAL|GUESS|SCF)\s*\(\S.+\)")
    

    __re_SCF_total_energy_section_begin = re.compile(">>>>\s*(\S*?Total\s*(Energy|energy)).*\[(\S.+)\]")

    # SCF section
    __re_SCF_start_section = re.compile(">>>>\s*(\S.+)\s+\[(\S.+)\]")
    __re_SCF_end_section = re.compile("<<<<\s+\[(\S.+)\]")

    # iteration
    __re_SCF_number_of_iteration = re.compile("(number of iteration|number_iteration)\s+=\s+(\d+)")

    # convcheck
    __re_SCF_convergence_check_section_begin = re.compile("(>>>>)?\s*\S*Conv.*(C|c)heck")

    # DfSummary section
    #__re_SCF_summary_section_begin = re.compile(">>>>\s*(\S*?Summary)\s+\[(\S.+)?\]")
    __re_SCF_summary_section_begin = re.compile("(>>>>)?\s*(Mulliken Population)")

    #
    # member function
    #
    def __init__(self, file_path = ""):
        """
        create the object

        Argument:
        file_path:      string of the ProteinDF output file path
        """
        self.__dfdata = DfData()

        self.__SCF = {} # dictionary type
        self.__SCF_current_section_name = ""
        self.__SCF_current_section_start_time = ""
        self.__SCF_current_iteration = 0
        self.__SCF_last_iteration = 0
        self.__stat = {}

        self.__step = ""

        self.__file_path = file_path
        self.__file = None # file object
        self.__read()


    def __str__(self):
        str = ""
        str += "num of atoms = %d\n" % (self.get_number_of_atoms())
        str += "num of AOs = %d\n" % (self.get_number_of_orbitals())
        str += "method = %s\n" % (self.get_method())
        if (self.get_method() == "RKS"):
            str += "occupation level(alpha) = %s\n" % (self.get_occupation_level("RKS"))
        elif (self.get_method() == "UKS"):
            str += "occupation level(alpha) = %s\n" % (self.get_occupation_level("UKS_ALPHA"))
            str += "occupation level(beta)  = %s\n" % (self.get_occupation_level("UKS_BETA"))
        str += "max iteration = %d\n" % (self.get_last_iteration())
        return str


    def get_dfdata(self):
        dfdata = copy.deepcopy(self.__dfdata)
        return dfdata
           
    #
    # private member function
    #
    def __read(self):
        """reading ProteinDF output file and parse."""
        self.__file = open(self.__file_path, "r")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break

            line = line.rstrip('\n')

            # step
            matchObj = self.__re_step.match(line);
            if (matchObj != None):
                self.__step = matchObj.group(1)
                self.__step = self.__step.upper()
                continue

            if (self.__re_input_data_start_section.search(line) != None):
                line = self.__read_input_data_section(line)

            # SCF step
            if (self.__step == "SCF"):
                # iteration
                self.__read_SCF_number_of_iteration(line)

                if (self.__re_SCF_total_energy_section_begin.search(line) != None):
                    line = self.__read_SCF_total_energy_section(line)

                if (self.__re_SCF_convergence_check_section_begin.search(line) != None):
                    line = self.__read_SCF_convcheck_section(line)
            
                if (self.__re_SCF_summary_section_begin.search(line) != None):
                    self.__read_SCF_summary_section(line)
                    continue

                # SCF section
                self.__read_SCF_start_section(line)
                self.__read_SCF_end_section(line)


        self.__file.close()


    def __read_input_data_section(self, line):
        """read INPUT DATA block"""
        re_input_keyword_value = re.compile("^\s*\[(\S.+)\]\s+\[\s*(\S.*?)\s*\]\s*$")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            if (re_input_keyword_value.search(line) != None):
                matchObj = re_input_keyword_value.match(line)
                keyword = matchObj.group(1)
                value = matchObj.group(2)

                if (keyword == 'control-number-of-atoms'):
                    self.__dfdata.set_number_of_atoms(value)
                elif (keyword == 'control-norb'):
                    self.__dfdata.set_number_of_orbitals(value)
                elif (keyword == 'method'):
                    self.__dfdata.set_method(value)
                elif (keyword == 'method/nsp/occlevel'):
                    self.__dfdata.set_occupation_level('RKS', value)
                elif (keyword == 'method/sp/alpha-spin-occlevel'):
                    self.__dfdata.set_occupation_level('UKS_ALPHA', value)
                elif (keyword == 'method/sp/beta-spin-occlevel'):
                    self.__dfdata.set_occupation_level('UKS_BETA', value)
                continue

            if (self.__re_input_data_end_section.search(line) != None):
                break

        return line


    def __read_SCF_total_energy_section(self, line):
        """read SCF total energy block"""
        self.__read_SCF_start_section(line)
        self._SCF_total_energy_section_name = self.__SCF_current_section_name
        re_SCF_total_energy_TE = re.compile("\s*TE\s*=\s*(\S+)")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            matchObj = re_SCF_total_energy_TE.match(line)
            if (matchObj != None):
                iteration = self.__SCF_current_iteration
                section_name = self.__SCF_current_section_name
                total_energy = matchObj.group(1)

                self.__dfdata.set_total_energy(
                    iteration = iteration,
                    total_energy = total_energy)

            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
                break

        line = self.__file.readline()
        return line

    def __read_SCF_convcheck_section(self, line):
        """read SCF convergence check block"""
        self.__read_SCF_start_section(line)
        re_SCF_convcheck_TE = re.compile(".*deviation\s*of\s*total\s*energy\s*=\s*(\S+)")
        re_SCF_convcheck_DensMat = re.compile(".*deviation\s*of\s*density\s*matrix\s*=\s*(\S+)")
        re_SCF_convcheck_CD = re.compile(".*deviation\s*of\s*cd\s*=\s*(\S+)")
        re_SCF_convcheck_KS = re.compile(".*deviation\s*of\s*kohn-sham\s*matrix\s*=\s*(\S+)")

        max_dev_TE = None
        max_dev_DensMat = None
        max_dev_CD = None
        max_dev_KS = None

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            matchObj = re_SCF_convcheck_TE.match(line)
            if (matchObj != None):
                max_dev_TE = float(matchObj.group(1))
                continue

            matchObj = re_SCF_convcheck_DensMat.match(line)
            if (matchObj != None):
                max_dev_DensMat = float(matchObj.group(1))
                continue

            matchObj = re_SCF_convcheck_CD.match(line)
            if (matchObj != None):
                max_dev_CD = float(matchObj.group(1))
                continue

            matchObj = re_SCF_convcheck_KS.match(line)
            if (matchObj != None):
                max_dev_KS = float(matchObj.group(1))
                continue

            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
                break

        # finalize
        iteration = self.__SCF_current_iteration
        self.__dfdata.set_convergence_info(
            iteration = iteration,
            max_deviation_of_total_energy = max_dev_TE,
            max_deviation_of_density_matrix = max_dev_DensMat,
            max_deviation_of_CD = max_dev_CD,
            max_deviation_of_KS_matrix = max_dev_KS
            )
        
        line = self.__file.readline()
        return line


    def __read_SCF_summary_section(self, line):
        """read SCF summary block"""
        self.__read_SCF_start_section(line)
        iteration = self.__SCF_current_iteration

        # variables
        ## ex)      Atom            tr(PS)         Net
        re_mulliken_atom_population_header = re.compile("^\s*(Atom|ATOM)\s+(Tr|tr)\(PS\)\s+Net");
        ## ex) N      1        7.690878      -0.690878
        re_mulliken_atom_population = re.compile("^\s*(\S+)\s+(\d+)\s+(\d+[\.]?\d*)\s+([+-]?\d+[\.]?\d*)")
        ## ex) GTO  ATOM  SHELL  tr(PS)
        re_mulliken_orbital_population_header = re.compile("^\s*GTO\s+(Atom|ATOM)\s+SHELL\s+(Tr|tr)\(PS\)")
        ## ex) 1   H      1   s()                    0.536552
        re_mulliken_orbital_population_start = re.compile("^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+([+-]?\d+[\.]?\d*)")
        ## ex) 2              s()                    0.087215
        re_mulliken_orbital_population_cont = re.compile("^\s*(\d+)\s+(\S+)\s+([+-]?\d+[\.]?\d*)")
        
        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            # mulliken atom population
            if (re_mulliken_atom_population_header.search(line) != None):
                # skip next blank line
                while (True):
                    line = self.__file.readline() # next line
                    line = line.rstrip('\n')
                    if (len(line) != 0):
                        break

                while (True):
                    if (re_mulliken_atom_population.search(line) != None):
                        matchObj = re_mulliken_atom_population.match(line)
                        symbol = matchObj.group(1)
                        atom_index = int(matchObj.group(2)) -1 # file shows the number begining from 1
                        gross_population = float(matchObj.group(3))
                        mulliken_population = float(matchObj.group(4))

                        self.__dfdata.set_mulliken_atom_population(
                            iteration = iteration,
                            atom_index = atom_index,
                            atom_symbol = symbol,
                            population = mulliken_population)
                    elif (len(line) != 0):
                        break

                    line = self.__file.readline()
                    line = line.rstrip('\n')

            # mulliken orbital population
            if (re_mulliken_orbital_population_header.search(line) != None):
                # skip next blank line
                while (True):
                    line = self.__file.readline() # next line
                    line = line.rstrip('\n')
                    if (len(line) != 0):
                        break

                gto_index = 0
                symbol = ""
                index = 0
                shell = ""
                gross_population = 0.0
                while (True):
                    if (re_mulliken_orbital_population_start.search(line) != None):
                        matchObj = re_mulliken_orbital_population_start.match(line)
                        gto_index = int(matchObj.group(1))
                        atom_symbol = matchObj.group(2)
                        atom_index = int(matchObj.group(3))
                        shell = matchObj.group(4)
                        gross_population = float(matchObj.group(5))
                    elif (re_mulliken_orbital_population_cont.search(line) != None):
                        matchObj = re_mulliken_orbital_population_cont.match(line)
                        gto_index = int(matchObj.group(1))
                        shell = matchObj.group(2)
                        gross_population = float(matchObj.group(3))
                    else:
                        check_blank_line = line.rstrip()
                        if (len(check_blank_line) == 0):
                            line = self.__file.readline()
                            line = line.rstrip('\n')
                            continue
                        else:
                            break
                

                    # set parameters
                    self.__dfdata.set_mulliken_orbital_population(
                        iteration = iteration,
                        gto_index = gto_index,
                        shell = shell,
                        population = gross_population
                        )

                    line = self.__file.readline()
                    line = line.rstrip('\n')

            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
                break
            
    def __read_SCF_start_section(self, line):
        matchObj = self.__re_SCF_start_section.match(line)
        if (matchObj != None):
            section_name = matchObj.group(1)
            section_start_time = matchObj.group(2)
            #print "SCF start section: %s %s" % (section_name, section_start_time)
            self.__SCF_current_section_name = section_name
            self.__SCF_current_section_start_time = section_start_time

    def __read_SCF_end_section(self, line):
        matchObj = self.__re_SCF_end_section.match(line)
        if (matchObj != None):
            section_end_time = matchObj.group(1)

            # prepare dictionary
            iteration = self.__SCF_current_iteration
            section_name = self.__SCF_current_section_name
            section_name = section_name.strip()

            start_time_str = self.__SCF_current_section_start_time
            end_time_str = section_end_time
            start_time = time.mktime(time.strptime(start_time_str, "%Y/%m/%d %H:%M:%S")) # 2007/08/20 22:42:52 format
            end_time = time.mktime(time.strptime(end_time_str, "%Y/%m/%d %H:%M:%S"))
            self.__dfdata.set_SCF_time(iteration, section_name, start_time, end_time)


    def __read_SCF_number_of_iteration(self, line):
        matchObj = self.__re_SCF_number_of_iteration.match(line)
        if (matchObj != None):
            iteration = matchObj.group(2) #matchObj.group(1)
            self.__SCF_current_iteration = int(iteration)
            self.__SCF_last_iteration = max(self.__SCF_last_iteration, self.__SCF_current_iteration)
            self.__dfdata.set_number_of_iterations(self.__SCF_last_iteration)

# end of class            
