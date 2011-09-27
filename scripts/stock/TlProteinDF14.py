#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import time
import xml.dom.minidom
import TlAtom

class PdfInput(object):
    # member variable ==================================================
    __file_path = ""
    __file = None # file object

    __current_group = ""
    __data = {}

    # comiled regulae explations =======================================
    __re_start_section = re.compile(">>>>\s*(\S+)\s*")
    __re_keyword_and_value = re.compile("\s*\[(.*)\]\s*\[(.*)\]")

    def get_groups(self):
        return self.__data.keys()

    def get_keywords(self, group):
        group = str(group)
        answer = None
        if (self.__data.has_key(group)):
            answer = self.__data[group].keys()
        return answer

    def get_value(self, group, keyword):
        group = str(group)
        keyword = str(keyword)
        answer = None
        if ((self.__data.has_key(group) == True) and
            (self.__data[group].has_key(keyword) == True)):
            answer = self.__data[group][keyword]
        return answer

    def __init__(self, file_path):
        self.__file_path = file_path
        self.__read()

    def __str__(self):
        output = ""
        for group in self.__data.keys():
            output += ">>>>%s\n" % (group)
            for keyword, value in self.__data[group].iteritems():
                output += "[%s] = [%s]\n" % (keyword, value)
        return output

    def __read(self):
        """read ProteinDF input file (e.g. fl_Input/fl_Globalinput) and parse it."""

        self.__file = open(self.__file_path, "r")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            if (self.__re_start_section.search(line) != None):
                self.__read_start_section(line)
                continue
            if (self.__re_keyword_and_value.search(line) != None):
                self.__read_keyword_and_value(line)
                continue

        self.__file.close()

    def __read_start_section(self, line):
        matchObj = self.__re_start_section.match(line)
        group = matchObj.group(1)
        self.__current_group = group

    def __read_keyword_and_value(self, line):
        matchObj = self.__re_keyword_and_value.match(line)
        keyword = str(matchObj.group(1))
        value = str(matchObj.group(2))
        group = str(self.__current_group)
        if (self.__data.has_key(group) != True):
            self.__data[group] = {}
        self.__data[group][keyword] = value

class PdfOutput(object):
    # member variable ==================================================
    __file_path = ""
    __file = None # file object

    __input_data = {}

    __SCF = {} # dictionary type
    __SCF_current_section_name = ""
    __SCF_current_section_start_time = ""
    __SCF_current_iteration = 0
    __SCF_last_iteration = 0

    __stat = {}

    #
    # compiled regular expressions =====================================
    # 
    
    # input section
    __re_input_data_start_section = re.compile(">>>>\s*Globalinput\s*Keyword List\s*<<<<")

    # integral section
    #  old style: ">>>>\s*INTEGRAL"
    __re_integral_start_section = re.compile("^\s*>>>>\s*INTEGRAL")

    __re_SCF_total_energy_section_begin = re.compile(">>>>\s*(\S*?Total(Energy|energy))(\s+\[(\S.+)\])?")

    # SCF section
    #  old style: ">>>>\s*(\S.+)\s*"
    __re_SCF_start_section = re.compile(">>>>\s*(\S.+)\s*")

    # end section
    #  old style: ">>>>\s*(\S.+)\s*"
    __re_SCF_end_section = re.compile(">>>>\s*(\S.+)\s*")

    __re_SCF_number_of_iteration = re.compile("number_iteration\s+=\s+(\d+)")

    # DfSummary section
    #  old style: ">>>>\s*(\S*?Summary)"
    __re_SCF_summary_section_begin = re.compile(">>>>\s*\S*?Summary")

    #
    # member function
    #
    def get_last_iteration(self):
        return self.__SCF_last_iteration

    def get_number_of_atoms(self):
        answer = None
        if (self.__input_data.has_key('control-number-of-atoms') == True):
            answer = int(self.__input_data['control-number-of-atoms'])
        return answer

    def get_number_of_orbitals(self):
        answer = None
        if (self.__input_data.has_key('control-norb') == True):
            answer = int(self.__input_data['control-norb'])
        return answer

    def get_total_energy(self, iteration):
        answer = None
        if (self.__SCF.has_key(iteration) == True):
            if (self.__SCF[iteration].has_key('TotalEnergy') == True):
                answer = self.__SCF[iteration]['TotalEnergy']
        return answer

    def get_mulliken_atom_population(self, atom_index, iteration =0):
        atom_index = int(atom_index)
        if (iteration == 0):
            iteration = self.get_last_iteration()
        iteration = int(iteration)

        atom_symbol = None
        gross_population = None
        mulliken_population = None

        if (self.__SCF[iteration].has_key('atom_population') == True):
            if (self.__SCF[iteration]['atom_population'].has_key(atom_index) == True):
                atom_symbol = self.__SCF[iteration]['atom_population'][atom_index]['atom_symbol']
                gross_population = float(self.__SCF[iteration]['atom_population'][atom_index]['gross_population'])
                mulliken_population = float(self.__SCF[iteration]['atom_population'][int(atom_index)]['mulliken_population'])

        return [atom_symbol, gross_population, mulliken_population]

    def get_mulliken_orbital_population(self, gto_index, iteration =0):
        gto_index = int(gto_index)
        if (iteration == 0):
            iteration = self.get_last_iteration()
        iteration = int(iteration)

        atom_index = None
        atom_symbol = None
        shell = None
        gross_population = None

        if (self.__SCF[iteration].has_key('orbital_population') == True):
            if (self.__SCF[iteration]['orbital_population'].has_key(gto_index) == True):
                atom_index = int(self.__SCF[iteration]['orbital_population'][gto_index]['atom_index'])
                atom_symbol = str(self.__SCF[iteration]['orbital_population'][gto_index]['atom_symbol'])
                shell = str(self.__SCF[iteration]['orbital_population'][gto_index]['shell'])
                gross_population = float(self.__SCF[iteration]['orbital_population'][gto_index]['gross_population'])

        return [atom_index, atom_symbol, shell, gross_population]

    def print_stat(self):
        self.__gather_stat()

        print "[SCF each iteration statics]"
        for iteration in range(1, (self.__SCF_last_iteration +1)):
            print "* #% 2d iteration:" % (iteration)
            for section_name in self.__stat['SCF'][iteration].keys():
                if ((self.__stat['SCF'][iteration][section_name].has_key('start_time') == True) and
                    (self.__stat['SCF'][iteration][section_name].has_key('end_time') == True)):
                    start_time = self.__stat['SCF'][iteration][section_name]['start_time']
                    end_time = self.__stat['SCF'][iteration][section_name]['end_time']
                    elapsed_time = self.__differ_time(start_time, end_time)

                    padding_space = " " * (25 - len(section_name))
                    print "  %s:%s%12.1f" % (section_name, padding_space, elapsed_time)
            print

        print "[SCF total statics]"
        print " for %2d iteration:" % (self.__SCF_last_iteration)
        for section_name in self.__stat['summary'].keys():
            elapsed_time =  self.__stat['summary'][section_name]
            average = elapsed_time / self.__SCF_last_iteration
            padding_space = " " * (25 - len(section_name))

            print "  %s:%sTotal=%12.1f,  Ave.=%12.1f" % (section_name, padding_space, elapsed_time, average)

    #
    # private member function
    #
    def __init__(self, file_path):
        self.__file_path = file_path
        self.__read()

    def __read(self):
        """reading ProteinDF output file and parse."""

        self.__file = open(self.__file_path, "r")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break

            line = line.rstrip('\n')

            if (self.__re_input_data_start_section.search(line) != None):
                line = self.__read_input_data_section(line)

            if (self.__re_SCF_total_energy_section_begin.search(line) != None):
                line = self.__read_SCF_total_energy_section(line)
            
            if (self.__re_SCF_summary_section_begin.search(line) != None):
                self.__read_SCF_summary_section(line)
                continue

            if (self.__re_SCF_start_section.search(line) != None):
                self.__read_SCF_start_section(line)
            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
            if (self.__re_SCF_number_of_iteration.search(line) != None):
                self.__read_SCF_number_of_iteration(line)

        self.__file.close()

    def __read_input_data_section(self, line):
        """read INPUT DATA block"""
        # old style: "^\s*\d+\s+(\S+)\s+(\S.+)$"
        # new style: "^\s*\[(\S.+)\]\s+\[(\S.+)\]"
        re_input_keyword_value_old = re.compile("^\s*\d+\s+(\S+)\s+(\S.+)\s*$")
        re_input_keyword_value_new = re.compile("^\s*\[(\S.+)\]\s+\[(\S.+)\]\s*$")
        re_oldtype_globalinput_section = re.compile("^\s*>>>>\s*Globalinput\s*Keyword\s*List\s*<<<<\s*$")

        while True:
            line = self.__file.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')

            # newtype
            if (re_input_keyword_value_new.search(line) != None):
                matchObj = re_input_keyword_value_new.match(line)
                keyword = matchObj.group(1)
                value = matchObj.group(2)
                self.__input_data[keyword] = value
                continue

            # oldtype
            if (re_input_keyword_value_old.search(line) != None):
                matchObj = re_input_keyword_value_old.match(line)
                keyword = matchObj.group(1)
                value = matchObj.group(2)
                self.__input_data[keyword] = value

            if (self.__re_integral_start_section.search(line) != None):
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

            if (re_SCF_total_energy_TE.search(line) != None):
                matchObj = re_SCF_total_energy_TE.match(line)
                iteration = self.__SCF_current_iteration
                section_name = self.__SCF_current_section_name
                total_energy = matchObj.group(1)

                if (self.__SCF.has_key(iteration) == False):
                    self.__SCF[iteration] = {}
                if (self.__SCF[iteration].has_key(section_name) == False):
                    self.__SCF[iteration][section_name] = {}
                self.__SCF[iteration]['TotalEnergy'] = float(total_energy)

            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
                break

        line = self.__file.readline()
        return line

    def __read_SCF_summary_section(self, line):
        """read SCF Summary block"""
        self.__read_SCF_start_section(line)
        iteration = self.__SCF_current_iteration

        # variables
        re_mulliken_atom_population_header = re.compile("^\s*ATOM\s+Tr\(PS\)\s+Net"); # ATOM  Tr(PS)  Net
        ## ex) N      1        7.690878      -0.690878
        re_mulliken_atom_population = re.compile("^\s*(\S+)\s+(\d+)\s+(\d+[\.]?\d*)\s+([+-]?\d+[\.]?\d*)")
        ## ex) GTO  ATOM  SHELL  Tr(PS)
        re_mulliken_orbital_population_header = re.compile("^\s*GTO\s+ATOM\s+SHELL\s+Tr\(PS\)")
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
                        atom_index = int(matchObj.group(2))
                        gross_population = float(matchObj.group(3))
                        mulliken_population = float(matchObj.group(4))
                        
                        if (self.__SCF.has_key(iteration) == False):
                            self.__SCF[iteration] = {}
                        if (self.__SCF[iteration].has_key('atom_population') == False):
                            self.__SCF[iteration]['atom_population'] = {}
                        if (self.__SCF[iteration]['atom_population'].has_key(atom_index) == False):
                            self.__SCF[iteration]['atom_population'][atom_index] = {}

                        self.__SCF[iteration]['atom_population'][atom_index]['atom_symbol'] = symbol
                        self.__SCF[iteration]['atom_population'][atom_index]['gross_population'] = float(gross_population)
                        self.__SCF[iteration]['atom_population'][atom_index]['mulliken_population'] = float(mulliken_population)
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
                    if (self.__SCF.has_key(iteration) == False):
                        self.__SCF[iteration] = {}
                    if (self.__SCF[iteration].has_key('orbital_population') == False):
                        self.__SCF[iteration]['orbital_population'] = {}
                    if (self.__SCF[iteration]['orbital_population'].has_key(gto_index) == False):
                        self.__SCF[iteration]['orbital_population'][gto_index] = {}

                    self.__SCF[iteration]['orbital_population'][gto_index]['atom_index'] = atom_index
                    self.__SCF[iteration]['orbital_population'][gto_index]['atom_symbol'] = atom_symbol
                    self.__SCF[iteration]['orbital_population'][gto_index]['shell'] = shell
                    self.__SCF[iteration]['orbital_population'][gto_index]['gross_population'] = gross_population

                    line = self.__file.readline()
                    line = line.rstrip('\n')

            if (self.__re_SCF_end_section.search(line) != None):
                self.__read_SCF_end_section(line)
                break
            
    def __read_SCF_start_section(self, line):
        matchObj = self.__re_SCF_start_section.match(line)
        if (matchObj != None):
            section_name = matchObj.group(2)
            section_start_time = matchObj.group(3)

            self.__SCF_current_section_name = section_name
            self.__SCF_current_section_start_time = section_start_time

    def __read_SCF_end_section(self, line):
        matchObj = self.__re_SCF_end_section.match(line)
        if (matchObj != None):
            section_end_time = matchObj.group(1)

            # prepare dictionary
            iteration = self.__SCF_current_iteration
            section_name = self.__SCF_current_section_name
            if (self.__stat.has_key('SCF') != True):
                self.__stat['SCF'] = {}
            if (self.__stat['SCF'].has_key(iteration) != True):
                self.__stat['SCF'][iteration] = {}
            if (self.__stat['SCF'][iteration].has_key(section_name) != True):
                self.__stat['SCF'][iteration][section_name] = {}

            self.__stat['SCF'][iteration][section_name]['start_time'] = self.__SCF_current_section_start_time
            self.__stat['SCF'][iteration][section_name]['end_time'] = section_end_time

    def __read_SCF_number_of_iteration(self, line):
        matchObj = self.__re_SCF_number_of_iteration.match(line)
        if (matchObj != None):
            iteration = matchObj.group(1)
            self.__SCF_current_iteration = int(iteration)
            self.__SCF_last_iteration = max(self.__SCF_last_iteration, self.__SCF_current_iteration)

    def __gather_stat(self):
        for iteration in range(1, (self.__SCF_last_iteration +1)):
            for section_name in self.__stat['SCF'][iteration].keys():
                if ((self.__stat['SCF'][iteration][section_name].has_key('start_time') == True) and
                    (self.__stat['SCF'][iteration][section_name].has_key('end_time') == True)):
                    start_time = self.__stat['SCF'][iteration][section_name]['start_time']
                    end_time = self.__stat['SCF'][iteration][section_name]['end_time']
                    elapsed_time = self.__differ_time(start_time, end_time)

                    if (self.__stat.has_key('summary') != True):
                        self.__stat['summary'] = {}
                    if (self.__stat['summary'].has_key(section_name) != True):
                        self.__stat['summary'][section_name] = float(0.0)
                    self.__stat['summary'][section_name] += elapsed_time
        #

    def __differ_time(self, start_time_str, end_time_str):
        start_time = time.mktime(time.strptime(start_time_str, "%Y/%m/%d %H:%M:%S")) # 2007/08/20 22:42:52 format
        end_time = time.mktime(time.strptime(end_time_str, "%Y/%m/%d %H:%M:%S"))
        
        time_diff = end_time - start_time
        return time_diff

# end of class            

# ======================================================================
# PDF xml
#
class PdfXml(object):
    # variables
    __dom = None
    __number_of_atoms = None
    __number_of_orbitals = None
    __SCF_last_iteration = None
    __SCF = {}

    #
    def set_data_from_pdf_output(self, pdf_output):
        assert (isinstance(pdf_output, PdfOutput) == True)
        self.__number_of_atoms = pdf_output.get_number_of_atoms()
        self.__number_of_orbitals = pdf_output.get_number_of_orbitals()
        # SCF
        self.__SCF_last_iteration = pdf_output.get_last_iteration()
        for itr in range(1, (self.__SCF_last_iteration +1)):
            if (self.__SCF.has_key(itr) != True):
                self.__SCF[itr] = {}
            # total energy
            if (pdf_output.get_total_energy(itr) != None):
                self.__SCF[itr]['total_energy'] = float(pdf_output.get_total_energy(itr))
            # atom pop
            if (self.__number_of_atoms != None):
                for atom_index in range(1, (self.__number_of_atoms +1)):
                    atom_index = int(atom_index)
                    atom_symbol, gross_population, mulliken_population = pdf_output.get_mulliken_atom_population(atom_index, itr)
                    if (mulliken_population != None):
                        if (self.__SCF[itr].has_key('atom_population') != True):
                            self.__SCF[itr]['atom_population'] = {}
                        self.__SCF[itr]['atom_population'][atom_index] = mulliken_population

    #def set_eigenvector(self, pdf_eigenvector):
        
    def set_number_of_atoms(self, orbitals):
        self.__number_of_atoms = orbitals
        
    def get_number_of_atoms(self):
        return self.__number_of_atoms

    def set_number_of_orbitals(self, orbitals):
        self.__number_of_orbitals = orbitals
        
    def get_number_of_orbitals(self):
        return self.__number_of_orbitals

    def get_number_of_iterations(self):
        return self.__SCF_last_iteration

    def get_total_energy(self, iteration):
        answer = None
        iteration = int(iteration)
        if (self.__SCF.has_key(iteration) == True):
            if (self.__SCF[iteration].has_key('total_energy') == True):
                answer = self.__SCF[iteration]['total_energy']
        return answer

    # hidden functions
    def __init__(self, xml_file_path=""):
        if (xml_file_path == ""):
            impl = xml.dom.minidom.getDOMImplementation()
            self.__dom = impl.createDocument(None, "pdf", None)
        else:
            self.__dom = xml.dom.minidom.parse(file(xml_file_path))
            self.__set_data_from_dom()

    def __del__(self):
        self.__dom.unlink()

    def __str__(self):
        text = ""
        if (self.__dom != None):
            self.build_xml()
            #text = self.__dom.toxml()
            text = self.__dom.toprettyxml()
            return text

    def __set_data_from_dom(self):
        pdf = self.__dom.getElementsByTagName("pdf")[0]
        # model
        model = pdf.getElementsByTagName("model")[0]
        dimensions = model.getElementsByTagName("dimensions")[0]
        if (dimensions.attributes.has_key('atoms') == True):
            self.__number_of_atoms = int(dimensions.attributes['atoms'].value)
            #print self.__number_of_atoms
        if (dimensions.attributes.has_key('orbitals') == True):
            self.__number_of_orbitals = int(dimensions.attributes['orbitals'].value)
            #print self.__number_of_orbitals
        # scf
        scf = pdf.getElementsByTagName("scf")[0]
        if (scf.attributes.has_key("iterations") == True):
            self.__SCF_last_iteration = int(scf.attributes["iterations"].value)
        iterations = scf.getElementsByTagName("iteration")
        index = 0
        for iteration in iterations:
            cycle = int(iteration.attributes["cycle"].value)
            total_energy = iteration.getElementsByTagName("total_energy")[0]
            if (total_energy != None):
                if (total_energy.attributes.has_key("value") == True):
                    if (self.__SCF.has_key(cycle) != True):
                        self.__SCF[cycle] = {}
                    self.__SCF[cycle]['total_energy'] = float(total_energy.attributes["value"].value)
                    #print "%d th = %f" % (cycle, self.__SCF[cycle]['total_energy'])

    def build_xml(self):
        pdf = self.__dom.documentElement
        # model
        model = self.__get_model_element()
        pdf.appendChild(model)

        # integral
        integral = self.__dom.createElement('integral')
        pdf.appendChild(integral)

        # guess
        guess = self.__dom.createElement('guess')
        pdf.appendChild(guess)

        # scf
        scf = self.__get_scf_element()
        pdf.appendChild(scf)

    def __get_model_element(self):
        """build and return model xml element"""
        model = self.__dom.createElement('model')
        # dimensions
        dimensions = self.__dom.createElement('dimensions')

        if (self.__number_of_atoms != None):
            dimensions.setAttribute('atoms', str(self.__number_of_atoms))
        if (self.__number_of_orbitals != None):
            dimensions.setAttribute('orbitals', str(self.__number_of_orbitals))

        model.appendChild(dimensions)
        return model

    def __get_scf_element(self):
        """build and return SCF xml element"""
        scf = self.__dom.createElement('scf')
        if (self.__SCF_last_iteration != None):
            scf.setAttribute('iterations', str(self.__SCF_last_iteration))

            # each iteration
            for itr in range(1, (self.__SCF_last_iteration +1)):
                if (self.__SCF.has_key(itr) == True):
                    itr_element = self.__dom.createElement('iteration')
                    itr_element.setAttribute('cycle', str(itr))

                    # total energy
                    if (self.__SCF[itr].has_key('total_energy') == True):
                        total_energy_element = self.__dom.createElement('total_energy')
                        total_energy = self.__SCF[itr]['total_energy']
                        total_energy_element.setAttribute('value', str(total_energy))
                        itr_element.appendChild(total_energy_element)

                    # population
                    if (self.__SCF[itr].has_key('atom_population') == True):
                        population_element = self.__dom.createElement('population')
                        population_element.setAttribute('type', 'atom')
                        for atom_index in self.__SCF[itr]['atom_population'].keys():
                            atom_index = int(atom_index)
                            atom_element = self.__dom.createElement('atom')
                            atom_element.setAttribute("index", str(atom_index))
                            value = self.__SCF[itr]['atom_population'][atom_index]
                            atom_element.setAttribute("value", str(value))
                            population_element.appendChild(atom_element)
                        itr_element.appendChild(population_element)
                        
                scf.appendChild(itr_element)

        return scf

    def get_xml(self):
        return self.__str__()

    def sysout(self):
        if (self.__dom != None):
            self.__dom.writexml(sys.stdout)

# end of class

class ProteinDFParser:
    # member variables
    _twoei = {}
    # comiled regulae explations =======================================
    #
    #  for Gross Atom Population
    #   ex)     N      1        7.690878      -0.690878
    _re_gross_atom_population = re.compile("\s*(\S+)\s+(\d+)\s+(\d+[\.]?\d*)\s+([+-]?\d+[\.]?\d*)")
    #  for J matrix
    _re_j_matrix = re.compile("J\(\s*(\d+),\s*(\d)\)\s*=\s*(\S+)")
    #  for TwoElectronIntegrals
    _re_two_electron_integrals = re.compile("I=\s*(\d+)\s*J=\s*(\d+)\s*K=\s*(\d+)\s*L=\s*(\d+)\s*Int=\s*([-]?\d\.\d+D[+|-]\d+)")
    
    # member functions =================================================
    def __init__(self, file_path):
        self.file_path = file_path

    def set_file_path(self, file_path):
        self.file_path = file_path

    def get_file_path(self):
        return self.file_path

    def get_twoei(self):
        return self._twoei
    
    # read and parse "fl_Out_Std"
    def read(self):
        re_GAP = re.compile("Gross atom population");
        fi = open(self.file_path, "r")
        for line in fi.readlines():
            line = line.rstrip('\n')

            if (re_GAP.search(line) != None):
                # read Gross Atom Population Block
                line = self.parse_gross_atom_population(fi)
        
            self._parse_two_electron_integrals(line)
        fi.close()

    # parse Gross Atom Population block.
    # if the end of block is found, return next line
    def parse_gross_atom_population(self, filein):
        for line in filein.readlines():
            line = line.rstrip('\n')
            if (line == ""):
                continue
            matchObj = self._re_gross_atom_population(line)
            if (matchObj != None):
                atom_symbol  = matchObj.group(1)
                atom_index = int(matchObj.group(2))
                atom_gross_population = float(matchObj.group(3))
                atom_mulliken_population = fload(matchObj.group(4))

                atom = TlAtom.TlAtom(atom_index, atom_symbol)
                atom.set_charge(atom_mulliken_population)
                
            else:
                break
        return line


    def _parse_two_electron_integrals(self, line):
        matchObj = self._re_two_electron_integrals.search(line)
        if (matchObj != None):
            i = int(matchObj.group(1))
            j = int(matchObj.group(2))
            k = int(matchObj.group(3))
            l = int(matchObj.group(4))
            t = (i, j, k, l)
            v = float(matchObj.group(5).replace("D", "E"))
            self._twoei[t] = v
            #print "(%3d %3d | %3d %3d) = % e" % (i, j, k, l, v)

    def _parse_J_matrix(self, line):
        matchObj = self._re_j_matrix.search(line)
        if (matchObj != None):
            x = int(matchObj.group(1))
            y = int(matchObj.group(2))
            v = float(matchObj.group(3))

    def _parse_matrix(self, file_in):
        re_comment = re.compile("^\s*[-]*$");
        for line in file_in.readlines():
            line = line.rstrip('\n')
            
            if (re_comment.search(line) != None):
                continue

            
            
