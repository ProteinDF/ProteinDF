#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import time
import xml.dom.minidom
import TlAtom
import TlXml

class PdfInput(object):
    """
    parse and analize ProteinDF input file (eg. fl_Globalinput)
    """

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




class PdfXml(object):
    """
    create and parse ProteinDF XML
    """

    # variables
    __dom = None
    __number_of_atoms = None
    __number_of_orbitals = None
    __SCF_last_iteration = None
    __SCF = {}

    def __init__(self, xml_file_path=""):
        """
        create the XML object from the XML file; otherwise empty object.
        
        Argument:
        xml_file_path:  string of the XML file path
        """
        if (xml_file_path == ""):
            impl = xml.dom.minidom.getDOMImplementation()
            self.__dom = impl.createDocument(None, "pdf", None)
        else:
            self.__dom = xml.dom.minidom.parse(file(xml_file_path))
            self.__set_data_from_dom()

    #
    def set_data_from_pdf_output(self, pdf_output):
        """
        set data from the PdfOutput object.

        Argument:
        pdf_output:     the object of PdfOutput.
        """
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
        
    def set_number_of_atoms(self, atoms):
        """set the number of atoms.
        
        Argument:
        atoms:          the number of atoms.
        """
        self.__number_of_atoms = int(atoms)
        
    def get_number_of_atoms(self):
        """return the number of atoms; otherwise returns None."""
        return self.__number_of_atoms

    def set_number_of_orbitals(self, orbitals):
        """set the number of orbitals.
        
        Argument:
        orbitals:       the number of orbitals.
        """
        self.__number_of_orbitals = orbitals
        
    def get_number_of_orbitals(self):
        """return the number of orbitals; otherwise returns None."""
        return self.__number_of_orbitals

    def get_number_of_iterations(self):
        """return the number of iterations; otherwise returns None."""
        return self.__SCF_last_iteration

    def get_total_energy(self, iteration):
        """
        return the total energy at the iteration; otherwise returns None.

        Argument:
        iteration:      the number of the iteration.
        """
        answer = None
        iteration = int(iteration)
        if (self.__SCF.has_key(iteration) == True):
            if (self.__SCF[iteration].has_key('total_energy') == True):
                answer = self.__SCF[iteration]['total_energy']
        return answer

    def get_mulliken_atom_population(self, atom_index, iteration):
        """
        return the Mulliken atom population at the iteration; otherwise returns None.

        Argument:
        atom_index:     the index of the atom.
        iteration:      the number of the iteration.
        """
        answer = None
        atom_index = int(atom_index)
        iteration = int(iteration)
        if (self.__SCF.has_key(iteration) == True):
            answer = 1
        return answer

    # hidden functions
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
            total_energy = iteration.getElementsByTagName('total_energy')[0]
            if (total_energy != None):
                if (total_energy.attributes.has_key("value") == True):
                    self.__SCF.setdefault(cycle, {})
                    self.__SCF[cycle]['total_energy'] = float(total_energy.attributes["value"].value)
            population = iteration.getElementsByTagName('population')
            if (len(population) != 0):
                atom_populations = population.getElementsByTagName("atom")
                self.__SCF.setdefault(cycle, {})
                self.__SCF[cycle].setdefault('atom_population', {})

                for atom_pop in atom_populations:
                    index = atom_pop.attributes.get('index', None)
                    value = atom_pop.attributes.get('value', None)
                    if ((index != None) and (value != None)):
                        index = int(index)
                        value = float(value)
                        self.__SCF[cycle]['atom_population'][index] = value
                        print "%d th SCF atom=%d pop=%f" % (cycle, index, value)
                

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
                        #population_element.setAttribute('type', 'atom')
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

            
            
