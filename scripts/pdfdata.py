#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import xml.dom.minidom

class PdfData(object):
    def __init__(self, data = "", options = {}):
        # variable initialize
        self.model = {}
        self.scf = {}
        self.logLevel = options.setdefault('log_level', 0)

    # model
    def set_number_of_atoms(self, atoms):
        atoms = int(atoms)
        self.model['number_of_atoms'] = atoms

    def get_number_of_atoms(self):
        return self.model.get('number_of_atoms', None)

    def set_number_of_orbitals(self, orbitals):
        orbitals = int(orbitals)
        self.model['number_of_orbitals'] = orbitals

    def get_number_of_orbitals(self):
        return self.model.get('number_of_orbitals', None)

    def set_method(self, method):
        method = str(method)
        self.model['method'] = method

    def get_method(self):
        return self.model.get('method', 'NSP')

    def set_occupation_level(self, spin_type, level):
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        self.model.setdefault('occupation_level', {})
        self.model['occupation_level'][spin_type] = level

    def get_occupation_level(self, spin_type):
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        level = None
        if (self.model.has_key('occupation_level') == True):
            if (self.model['occupation_level'].has_key(spin_type) == True):
                level = self.model['occupation_level'][spin_type]
        return level


    def get_HOMO_level(self, method):
        level = self.get_occupation_level(method)
        matchObj = re.match('(\d+)\s*-\s*(\d+)', level)
        homo_level = None
        if (matchObj != None):
            homo_level = int(matchObj.group(2))
        return homo_level


    # SCF
    def set_number_of_iterations(self, iterations):
        self.scf['number_of_iterations'] = iterations

    def get_number_of_iterations(self):
        return self.scf.get('number_of_iterations', None)

    def set_total_energy(self, iteration, total_energy):
        iteration = int(iteration)
        self.scf.setdefault(iteration, {})
        self.scf[iteration]['total_energy'] = total_energy

    def get_total_energy(self, iteration = 0):
        if (iteration == 0):
            iteration = self.get_number_of_iterations()
        return self.scf[int(iteration)].get('total_energy', None)

    def set_mulliken_atom_population(self, iteration, atom_index, population, atom_symbol =""):
        iteration = int(iteration)
        atom_index = int(atom_index)
        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('mulliken_atom_population', {})
        self.scf[iteration]['mulliken_atom_population'].setdefault(atom_index, {})
        self.scf[iteration]['mulliken_atom_population'][atom_index]['atom_symbol'] = atom_symbol
        self.scf[iteration]['mulliken_atom_population'][atom_index]['population'] = population

    def get_mulliken_atom_population(self, iteration, atom_index):
        iteration = int(iteration)
        atom_index = int(atom_index)
        population = None
        if (self.scf.has_key(iteration) == True):
            if (self.scf[iteration].has_key('mulliken_atom_population') == True):
                mulliken_pop = self.scf[iteration]['mulliken_atom_population']
                if (mulliken_pop.has_key(atom_index) == True):
                    population = mulliken_pop[atom_index].get('population')
        return population


    def set_mulliken_orbital_population(self, iteration, gto_index, shell, population):
        iteration = int(iteration)
        gto_index = int(gto_index)
        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('mulliken_orbital_population', {})
        self.scf[iteration]['mulliken_orbital_population'].setdefault(gto_index, {})
        self.scf[iteration]['mulliken_orbital_population'][gto_index]['shell'] = shell
        self.scf[iteration]['mulliken_orbital_population'][gto_index]['population'] = population


    def get_mulliken_orbital_population(self, iteration, gto_index):
        iteration = int(iteration)
        gto_index = int(gto_index)
        population = None
        if (self.scf.has_key(iteration) == True):
            if (self.scf[iteration].has_key('mulliken_orbital_population') == True):
                mulliken_pop = self.scf[iteration]['mulliken_orbital_population']
                if (mulliken_pop.has_key(gto_index) == True):
                    population = mulliken_pop[gto_index].get('population')
        return population


    def set_energy_levels(self, iteration, spin_type, energy_levels, unit = 'a.u.'):
        '''
        set energy level data by iteration step.
        energy_levels needs to separate white space.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert(unit == 'a.u.')

        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('energy_levels', {})
        self.scf[iteration]['energy_levels'][spin_type] = str(energy_levels)
        self.scf[iteration]['energy_levels']['unit'] = unit


    def get_energy_levels(self, iteration = 0, spin_type = 'RKS', unit = 'a.u.'):
        '''
        get energy level data by iteration step.
        energy_levels separate white space.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert(unit == 'a.u.')

        answer = None
        if (iteration == 0):
            iteration = self.get_number_of_iterations()
        energy_levels = self.scf[iteration].get('energy_levels', None)
        if (energy_levels != None):
            data = self.scf[iteration]['energy_levels'].get(spin_type, None)
            if (data != None):
                answer = data
        return answer


    def get_band_gap(self, iteration = 0, spin_type = 'RKS', unit = 'a.u.'):
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert(unit == 'a.u.')
        
        energy_levels = self.get_energy_levels(iteration, spin_type, unit)
        HOMO_level = self.get_HOMO_level(spin_type)
        energy_level_list = energy_levels.split(",")
        # 1au = 27.21138 eV
        HOMO_energy = energy_level_list[HOMO_level -1]; # '-1' is needed by option base 0.
        LUMO_energy = energy_level_list[HOMO_level]; # '+1 -1' = 0
        return float(LUMO_energy - HOMO_energy)


    def __str__(self):
        out = ""
        if (self.get_number_of_atoms() != None):
            out += "# of atoms: %d\n" % (self.get_number_of_atoms())
        if (self.get_number_of_orbitals() != None):
            out += "# of orbitals: %d\n" % (self.get_number_of_orbitals())

        method = self.get_method()
        if (method != None):
            out += "method: %s\n" % (method)
            if (method.upper() == "NSP"):
                occ = self.get_occupation_level(spin_type = 'RKS')
                if (occ != None):
                    out += "occupation level: %s\n" % (occ)
            elif (method.upper() == "SP"):
                occ_a = self.get_occupation_level(spin_type = 'UKS_ALPHA')
                if (occ_a != None):
                    out += "occupation level(alpha): %s\n" % (occ_a)
                occ_b = self.get_occupation_level(spin_type = 'UKS_BETA')
                if (occ_b != None):
                    out += "occupation level(beta): %s\n" % (occ_b)

        iterations = self.get_number_of_iterations()
        for itr in range(1, (iterations +1)):
            out += "iteration step = %d:\n" % (itr)
            total_energy = self.get_total_energy(itr)
            if (total_energy != None):
                out += " total energy = %f\n" % (float(total_energy))
            method = self.get_method()
            if (method.upper() == "NSP"):
                e_levels = self.get_energy_levels(
                    iteration = itr,
                    spin_type = 'RKS')
                if (e_levels != None):
                    out += " energy levels: %s\n" % (e_levels)
            elif (method.upper() == "SP"):
                e_levels_a = self.get_energy_levels(
                    iteration = itr,
                    spin_type = 'UKS_ALPHA')
                if (e_levels_a != None):
                    out += " energy levels(alpha): %s\n" % (e_levels_a)
                e_levels_b = self.get_energy_levels(
                    iteration = itr,
                    spin_type = 'UKS_BETA')
                if (e_levels_b != None):
                    out += " energy levels(alpha): %s\n" % (e_levels_b)

        return out


    def setByXmlFile(self, xmlFilePath):
        dom = xml.dom.minidom.parse(file(xmlFilePath))
        self.setByXmlDom(dom)
        dom.unlink()

    def setByXmlString(self, xmlString):
        dom = xml.dom.minidom.parseString(xmlString)
        self.setByXmlDom(dom)
        dom.unlink()

    def setByXmlDom(self, dom):
        pdf = dom.getElementsByTagName("pdf")[0]
        # model
        model = pdf.getElementsByTagName("model")[0]
        dimensions = model.getElementsByTagName("dimensions")[0]
        if (dimensions.attributes.has_key('atoms') == True):
            self.model['numberOfAtoms'] = int(dimensions.attributes['atoms'].value)
        if (dimensions.attributes.has_key('orbitals') == True):
            self.model['numberOfOrbitals'] = int(dimensions.attributes['orbitals'].value)
        # scf
        scf = pdf.getElementsByTagName("scf")[0]
        if (scf.attributes.has_key("iterations") == True):
            self.scf['numberOfIterations'] = int(scf.attributes["iterations"].value)
        iterations = scf.getElementsByTagName("iteration")
        for iteration in iterations:
            cycle = int(iteration.attributes["cycle"].value)
            self.scf.setdefault(cycle, {})

            if (len(iteration.getElementsByTagName('total_energy')) != 0):
                totalEnergy = iteration.getElementsByTagName('total_energy')[0]
                if (totalEnergy != None):
                    if (totalEnergy.attributes.has_key("value") == True):
                        self.scf[cycle]['totalEnergy'] = float(totalEnergy.attributes["value"].value)

            if (len(iteration.getElementsByTagName('energy_level')) != 0):
                for energyLevel in (iteration.getElementsByTagName('energy_level')):
                    spin_type = 'RKS'
                    if (energyLevel.attributes.has_key('type') == True):
                        spin_type = str(energyLevel.attributes['type'].value)
                        spin_type = spin_type.upper()
                    if (len(energyLevel.childNodes) != 0):
                        energyLevelData = energyLevel.childNodes[0]
                        self.scf[cycle].setdefault('energyLevel', {})
                        self.scf[cycle]['energyLevel'][spin_type] = str(energyLevelData.data)
            
            for population in iteration.getElementsByTagName('population'):
                atom_populations = population.getElementsByTagName("atom")
                self.scf.setdefault(cycle, {})
                self.scf[cycle].setdefault('mullikenAtomPopulation', {})

                for atom_pop in atom_populations:
                    index = atom_pop.attributes.get('index', None)
                    pop = atom_pop.attributes.get('value', None)
                    symbol = atom_pop.attributes.get('symbol', None)
                    if ((index != None) and (pop != None)):
                        i = int(index.value)
                        self.scf[cycle]['mullikenAtomPopulation'][i] = {'value': float(pop.value)}
                        if (symbol != None):
                            self.scf[cycle]['mullikenAtomPopulation'][i]['symbol'] = str(symbol.value)

# end of class

