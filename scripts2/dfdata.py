#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import array
import copy

from dfmatrix import *
from bratom import *

class DfData(object):
    def __init__(self, data = "", options = {}):
        self.model = {}
        self.scf = {}

    def update(self, data):
        """
        combine two DfData object.
        """
        if (isinstance(data, DfData) == True):
            self.model.update(data.model)
            self.scf.update(data.scf)

    # model
    def set_number_of_atoms(self, atoms):
        atoms = int(atoms)
        self.model['number_of_atoms'] = atoms


    def get_number_of_atoms(self):
        return self.model.get('number_of_atoms', None)


    def set_atom(self, index, atom):
        self.model.setdefault('coord', [])
        if (len(self.model['coord']) <= index):
            self.model['coord'].extend([None] * (index - len(self.model['coord']) +1))
        self.model['coord'][index] = atom.get_raw_data()


    def get_atom(self, index):
        atom = BrAtom()
        atom.set_by_raw_data(self.model['coord'][index])
        return atom

        
    def set_number_of_orbitals(self, orbitals):
        orbitals = int(orbitals)
        self.model['number_of_orbitals'] = orbitals


    def get_number_of_orbitals(self):
        return self.model.get('number_of_orbitals', None)


    def set_method(self, method):
        method = str(method)
        method = method.upper()
        # transcode
        if (method == 'NSP'):
            method = 'RKS'
        elif (method == 'SP'):
            method = 'UKS'
        assert((method == 'RKS') or (method == 'UKS') or (method == 'ROKS'))
        self.model['method'] = method


    def get_method(self):
        return self.model.get('method', 'RKS')


    def set_occupation_level(self, spin_type, level):
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        self.model.setdefault('occupation_level', {})
        level = level.strip()
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
        if (level == None):
            return None

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
        self.scf[iteration]['total_energy'] = float(total_energy)


    def get_total_energy(self, iteration = 0):
        if (iteration == 0):
            iteration = self.get_number_of_iterations()
        answer = self.scf[int(iteration)].get('total_energy', None)
        if (answer != None):
            answer = float(answer)
        return answer 


    def set_convergence_info(self, iteration = 0,
                             max_deviation_of_total_energy = None,
                             max_deviation_of_density_matrix = None,
                             max_deviation_of_CD = None,
                             max_deviation_of_KS_matrix = None
                             ):
        iteration = int(iteration)
        if (max_deviation_of_total_energy != None):
            self.scf.setdefault(iteration, {})
            self.scf[iteration].setdefault('convergence_info', {})
            self.scf[iteration]['convergence_info']['max_deviation_of_total_energy'] = float(max_deviation_of_total_energy)
        if (max_deviation_of_density_matrix != None):
            self.scf.setdefault(iteration, {})
            self.scf[iteration].setdefault('convergence_info', {})
            self.scf[iteration]['convergence_info']['max_deviation_of_density_matrix'] = float(max_deviation_of_density_matrix)
        if (max_deviation_of_CD != None):
            self.scf.setdefault(iteration, {})
            self.scf[iteration].setdefault('convergence_info', {})
            self.scf[iteration]['convergence_info']['max_deviation_of_CD'] = float(max_deviation_of_CD)
        if (max_deviation_of_KS_matrix != None):
            self.scf.setdefault(iteration, {})
            self.scf[iteration].setdefault('convergence_info', {})
            self.scf[iteration]['convergence_info']['max_deviation_of_KS_matrix'] = float(max_deviation_of_KS_matrix)


    def get_convergence_info(self, iteration = 0):
        if (iteration == 0):
            iteration = self.get_number_of_iteration()
        answer = {}
        if (self.scf.has_key(iteration) == True):
            if (self.scf[iteration].has_key('convergence_info') == True):
                for item in ('max_deviation_of_total_energy',
                             'max_deviation_of_density_matrix',
                             'max_deviation_of_CD',
                             'max_deviation_of_KS_matrix'
                             ):
                    if (self.scf[iteration]['convergence_info'].has_key(item) == True):
                        answer[item] = self.scf[iteration]['convergence_info'][item]

        return answer
            

    def set_mulliken_atom_population(self, iteration, atom_index, population):
        iteration = int(iteration)
        atom_index = int(atom_index)
        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('mulliken_atom_population', array.array('d'))
        mulliken_pop = self.scf[iteration]['mulliken_atom_population']
        prev_size = len(mulliken_pop)
        if (len(mulliken_pop) <= atom_index):
            mulliken_pop.extend([0.0] * (atom_index - len(mulliken_pop) +1))
        mulliken_pop[atom_index] = population

    def get_mulliken_atom_population(self, iteration, atom_index):
        iteration = int(iteration)
        atom_index = int(atom_index)
        population = None
        if (self.scf.has_key(iteration) == True):
            if (self.scf[iteration].has_key('mulliken_atom_population') == True):
                mulliken_pop = self.scf[iteration]['mulliken_atom_population']
                if (atom_index < len(mulliken_pop)):
                    population = mulliken_pop[atom_index]
                else:
                    print("get_mulliken_atom_population: out of atom_index")
            else:
                print("get_mulliken_atom_population: data not found ")
        else:
            print("get_mulliken_atom_population: iteration not found ")
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


    def set_lcao_matrix(self, iteration, spin_type, dfmatrix):
        '''
        set LCAO matrix data by iteration step.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert((isinstance(dfmatrix, DfMatrix) == True) or
               (isinstance(dfmatrix, DfSymmetricMatrix) == True))

        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('lcao', {})
        self.scf[iteration]['lcao'][spin_type] = dfmatrix.get_raw_data()


    def get_lcao_matrix(self, iteration, spin_type):
        '''
        get LACO matrix by iteration step.
        return the DfMatrix(DfSymmetricMatrix) object or None.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))

        answer = None
        if (iteration == 0):
            iteration = self.get_number_of_iterations()
        if (self.scf.has_key(iteration) == True):
            lcao = self.scf[iteration].get('lcao', None)
            if (lcao != None):
                raw_data = lcao.get(spin_type, None)
                if (raw_data != None):
                    matrix_type = raw_data.get('type', None)
                    if (matrix_type == 'RLHD'):
                        answer = DfSymmetricMatrix(raw_data)
                    elif (matrix_type == ''):
                        answer = DfMatrix(raw_data)

        return answer


    def set_energy_levels(self, iteration, spin_type, energy_levels, unit = 'a.u.'):
        '''
        set energy level data by iteration step.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert((isinstance(energy_levels, list) == True) or
               (isinstance(energy_levels, array.array) == True))
        assert(unit == 'a.u.')

        if (isinstance(energy_levels, array.array) == False):
            tmp = array('d')
            tmp.fromlist(energy_levels)
            energy_levels = tmp

        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault('energy_levels', {})
        self.scf[iteration]['energy_levels'][spin_type] = energy_levels
        self.scf[iteration]['energy_levels']['unit'] = unit


    def get_energy_levels(self, iteration = 0, spin_type = 'RKS', request_unit = 'a.u.'):
        '''
        get energy level data by iteration step.
        return the array object or None.
        '''
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))

        answer = None
        if (iteration == 0):
            iteration = self.get_number_of_iterations()
        if (self.scf.has_key(iteration) == True):
            energy_levels = self.scf[iteration].get('energy_levels', None)
            if (energy_levels == None):
                return None
            data = energy_levels.get(spin_type, None)
            if (data == None):
                return None
            #data = set(data)
            answer = copy.deepcopy(data)

            # data translate
            data_unit = self.scf[iteration]['energy_levels'].get('unit', 'a.u.')
            coef = 1.0 # a.u.を変換係数1.0とする
            if (data_unit != request_unit):
                if (data_unit == "eV"):
                    coef /= 27.21138
                if (request_unit == "eV"):
                    coef *= 27.21138

            tmp = []
            for level, value in enumerate(answer):
                #answer[level] = value * coef
                tmp.append(value * coef)
            answer = tmp
        
        return answer


    def get_band_gap(self, iteration = 0, spin_type = 'RKS', unit = 'a.u.'):
        iteration = int(iteration)
        spin_type = str(spin_type).upper()
        assert((spin_type == 'RKS') or (spin_type == 'UKS_ALPHA') or (spin_type == 'UKS_BETA'))
        assert(unit == 'a.u.')
        
        energy_levels = self.get_energy_levels(iteration, spin_type, unit)
        HOMO_level = self.get_HOMO_level(spin_type)
        #energy_level_list = energy_levels.split(",")
        # 1au = 27.21138 eV
        answer = None
        if (HOMO_level < len(energy_levels)):
            HOMO_energy = energy_levels[HOMO_level -1]; # '-1' is needed by option base 0.
            LUMO_energy = energy_levels[HOMO_level]; # '+1 -1' = 0
            answer = float(LUMO_energy - HOMO_energy)
        return answer

    def get_raw_data(self):
        raw = {}
        raw['model'] = self.model
        raw['scf'] = self.scf
        return raw

    def set_raw_data(self, raw):
        if (isinstance(raw, dict) == True):
            self.model = raw.setdefault('model', {})
            self.scf = raw.setdefault('scf', {})


    # ==========================================================================
    # for statics
    #
    def set_SCF_time(self, iteration, section_name, start_time, end_time):
        """
        """
        iteration = int(iteration)
        assert(isinstance(section_name, str) == True)
        section = section_name.upper()
        self.scf.setdefault(iteration, {})
        self.scf[iteration].setdefault("time", {})
        self.scf[iteration]["time"][section_name] = {"start_time": start_time,
                                                     "end_time": end_time}

    def get_SCF_time(self, iteration, section_name):
        answer = None
        if (self.scf.has_key(iteration) == True):
            if (self.scf[iteration].has_key("time") == True):
                if (self.scf[iteration]["time"].has_key(section_name) == True):
                    answer = [self.scf[iteration]["time"][section_name].get("start_time", 0.0),
                              self.scf[iteration]["time"][section_name].get("end_time", 0.0)]
        return answer
    
    
    def stat(self):
        answer = ""
        summary = {}
        total_time = 0.0
        number_of_iterations = self.get_number_of_iterations()

        answer += "[SCF each iteration statics]\n"
        for iteration in range(1, (number_of_iterations +1)):
            if ((self.scf.has_key(iteration) == False) or
                (self.scf[iteration].has_key("time") == False)):
                continue

            answer += "iteration: %d\n" % (iteration)
            for section_name in self.scf[iteration]["time"].keys():
                if (summary.has_key(section_name) == False):
                    summary.setdefault(section_name, 0.0)
                if ((self.scf[iteration]["time"][section_name].has_key('start_time') == True) and
                    (self.scf[iteration]["time"][section_name].has_key('end_time') == True)):
                    time_range = self.get_SCF_time(iteration, section_name)
                    start_time = time_range[0]
                    end_time = time_range[1]
                    delta_time = end_time - start_time
                    summary[section_name] += delta_time
                    total_time += delta_time
                    
                    padding_space = " " * (40 - len(section_name))
                    answer += "  %s:%s%12.1f\n" % (section_name, padding_space, delta_time)
            answer += "\n"

        answer += "[SCF total statics]"
        answer += " for %2d iteration:\n" % (number_of_iterations)
        for section_name, value in summary.iteritems():
            average = value / number_of_iterations
            in_percent = value / total_time * 100.0
            padding_space = " " * (40 - len(section_name))
            answer += "  %s:%sTotal=%12.1f,  Ave.=%12.1f (%5.2f%%)\n" % (
                section_name, padding_space, value, average, in_percent)

        return answer
    
        
    # ==========================================================================
    # string
    #
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

# end of class

