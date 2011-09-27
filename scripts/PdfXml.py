#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import re
import xml.dom.minidom

from pdfdata import PdfData

class PdfXml(object):
    """
    create and parse ProteinDF XML
    """
    # variables
    #__dom = None
    #__pdfdata = PdfData()
    #__number_of_atoms = None
    #__number_of_orbitals = None
    #__SCF_last_iteration = None
    #__SCF = {}

    def __init__(self, input = None):
        """
        create the XML object from a XML file or PdfData object
        if input is nothing, create empty object.
        
        Argument:
        input:  string of the XML file path, or PdfData object.
        """
        self.__initialize()
        if (isinstance(input, basestring) == True):
            self.__dom = xml.dom.minidom.parse(file(input))
            self.__set_data_from_dom()
        else:
            # build empty object
            impl = xml.dom.minidom.getDOMImplementation()
            self.__dom = impl.createDocument(None, "pdf", None)
            
            if (isinstance(input, PdfData) == True):
                self.set_data_from_pdfdata(input)

    def __initialize(self):
        self.__dom = None
        self.__pdfdata = PdfData()
        self.__number_of_atoms = None
        self.__number_of_orbitals = None
        self.__SCF_last_iteration = None
        self.__SCF = {}

    def set_data_from_pdfdata(self, pdfdata):
        """
        set data from the PdfOutput object.

        Argument:
        pdf_output:     object of PdfData.
        """
        assert (isinstance(pdfdata, PdfData) == True)
        self.__pdfdata = copy.deepcopy(pdfdata)


    def get_pdfdata(self):
        answer = copy.deepcopy(self.__pdfdata)
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
            atoms = int(dimensions.attributes['atoms'].value)
            self.__pdfdata.set_number_of_atoms(atoms)
        if (dimensions.attributes.has_key('orbitals') == True):
            orbitals = int(dimensions.attributes['orbitals'].value)
            self.__pdfdata.set_number_of_orbitals(orbitals)
            
        if (len(model.getElementsByTagName("method")) > 0):
            method = model.getElementsByTagName("method")[0]
            method_txt = self.__xml_get_text(method)
            method_txt = re.sub('^\s+', '', method_txt)
            method_txt = re.sub('\s*$', '', method_txt)
            self.__pdfdata.set_method(method_txt)
            
        occ_elements = model.getElementsByTagName("occupation_level")
        for occ in occ_elements:
            if (occ.attributes.has_key('type') == True):
                occ_spintype = occ.attributes['type'].value.upper()
                if (occ.attributes.has_key('value') == True):
                    occ_level = occ.attributes['value'].value
                    self.__pdfdata.set_occupation_level(
                        spin_type = occ_spintype,
                        level = occ_level)
                
        # scf
        scf = pdf.getElementsByTagName("scf")[0]
        if (scf.attributes.has_key("iterations") == True):
            iterations = int(scf.attributes["iterations"].value)
            self.__pdfdata.set_number_of_iterations(iterations)
        itr_elements = scf.getElementsByTagName("iteration")
        index = 0
        for iteration in itr_elements:
            step_obj = None
            if (iteration.attributes.has_key("step") == True):
                step_obj = iteration.attributes["step"]
            elif (iteration.attributes.has_key("cycle") == True):
                # for old version
                step_obj = iteration.attributes["cycle"]
            step = int(step_obj.value)

            # total energy
            if ((iteration.getElementsByTagName('total_energy') != None) and
                (len(iteration.getElementsByTagName('total_energy')) > 0)):
                total_energy = iteration.getElementsByTagName('total_energy')[0]
                if (total_energy != None):
                    if (total_energy.attributes.has_key("value") == True):
                        value = float(total_energy.attributes["value"].value)
                        self.__pdfdata.set_total_energy(
                            iteration = step, total_energy = value)

            # energy levels
            energy_levels = iteration.getElementsByTagName('energy_levels')
            if (energy_levels != None):
                for e_level in energy_levels:
                    if (e_level.attributes.has_key('type') == True):
                        spin = e_level.attributes['type'].value
                    else:
                        spin = 'RKS'
                    el = self.__xml_get_text(e_level)
                    el = re.sub('^\s+', '', el)
                    el = re.sub('\s*$', '', el)
                    self.__pdfdata.set_energy_levels(
                        iteration = step,
                        spin_type = spin,
                        energy_levels = el)

            # population
            if (len(iteration.getElementsByTagName('population')) != 0):
                population = iteration.getElementsByTagName('population')[0]
                if (population.nodeType == population.ELEMENT_NODE):
                    atom_populations = population.getElementsByTagName("atom")
                    if (atom_populations != None):
                        for atom_pop in atom_populations:
                            index = atom_pop.attributes.get('index', None)
                            value = atom_pop.attributes.get('value', None)
                            if ((index != None) and (value != None)):
                                index = int(index.value)
                                value = float(value.value)
                                self.__pdfdata.set_mulliken_atom_population(
                                    iteration = step,
                                    atom_index = index,
                                    population = value)

                
    def __xml_get_text(self, element):
        nodelist = element.childNodes
        rc = ""
        for node in nodelist:
            if (node.nodeType == node.TEXT_NODE or
                node.nodeType == node.CDATA_SECTION_NODE):
                rc = rc + str(node.data)
        return rc


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
        atoms = self.__pdfdata.get_number_of_atoms()
        if (atoms != None):
            dimensions.setAttribute('atoms', str(atoms))
        orbitals = self.__pdfdata.get_number_of_orbitals()
        if (orbitals != None):
            dimensions.setAttribute('orbitals', str(orbitals))
        model.appendChild(dimensions)

        # method
        method_el = self.__dom.createElement('method')
        method = self.__pdfdata.get_method()
        if (method != None):
            method_txt = self.__dom.createTextNode(method)
            method_el.appendChild(method_txt)
            model.appendChild(method_el)

        # occupation level
        method = self.__pdfdata.get_method()
        if (method == None):
            method = 'NSP'
        if (method.upper() == 'NSP'):
            occ_level = self.__pdfdata.get_occupation_level('RKS')
            if (occ_level != None):
                occ_el = self.__dom.createElement('occupation_level')
                occ_el.setAttribute('type', 'RKS')
                occ_el.setAttribute('value', occ_level)
                model.appendChild(occ_el)
        elif (method.upper() == 'SP'):
            occ_level_a = self.__pdfdata.get_occupation_level('UKS_ALPHA')
            if (occ_level_a != None):
                occ_el_a = self.__dom.createElement('occupation_level')
                occ_el_a.setAttribute('type', 'UKS_ALPHA')
                occ_el_a.setAttribute('value', occ_level_a)
                model.appendChild(occ_el_a)
            occ_level_b = self.__pdfdata.get_occupation_level('UKS_BETA')
            if (occ_level_b != None):
                occ_el_b = self.__dom.createElement('occupation_level')
                occ_el_b.setAttribute('type', 'UKS_BETA')
                occ_el_b.setAttribute('value', occ_level_b)
                model.appendChild(occ_el_b)

        return model


    def __get_scf_element(self):
        """build and return SCF xml element"""
        scf = self.__dom.createElement('scf')

        iterations = self.__pdfdata.get_number_of_iterations()
        if (iterations != None):
            scf.setAttribute('iterations', str(iterations))

            method = self.__pdfdata.get_method()
            if (method.upper() == 'NSP'):
                spin_list = ['RKS']
            elif (method.upper() == 'SP'):
                spin_list = ['UKS_ALPHA', 'UKS_BETA']

            # each iteration
            for itr in range(1, (iterations +1)):
                itr_element = self.__dom.createElement('iteration')
                itr_element.setAttribute('step', str(itr))

                # total energy
                total_energy = self.__pdfdata.get_total_energy(iteration = itr)
                if (total_energy != None):
                    total_energy_element = self.__dom.createElement('total_energy')
                    total_energy_element.setAttribute('value', str(total_energy))
                    itr_element.appendChild(total_energy_element)

                # energy level
                for spin in spin_list:
                    energy_level = self.__pdfdata.get_energy_levels(
                        iteration = itr, spin_type = spin)
                    if (energy_level != None):
                        el_element = self.__dom.createElement('energy_levels')
                        el_element.setAttribute('type', str(spin))
                        el_body = self.__dom.createCDATASection(str(energy_level));
                        el_element.appendChild(el_body);
                        itr_element.appendChild(el_element)
                
                # population
                atoms = self.__pdfdata.get_number_of_atoms()
                if (atoms != None):
                    population_element = self.__dom.createElement('population')
                    for index in range(0, atoms):
                        population = self.__pdfdata.get_mulliken_atom_population(
                            iteration = itr, atom_index = index)
                        if (population != None): 
                            atom_element = self.__dom.createElement('atom')
                            atom_element.setAttribute("index", str(index))
                            atom_element.setAttribute("value", str(population))
                            population_element.appendChild(atom_element)
                    itr_element.appendChild(population_element)

                # append iteration element
                scf.appendChild(itr_element)

        return scf


    def get_xml(self):
        return self.__str__()


    def sysout(self):
        if (self.__dom != None):
            self.__dom.writexml(sys.stdout)

# end of class

