import os
import re
import pandas as pd
from pyBigstick.orbits import df_orbits
import numpy as np

nuclear_data = pd.read_json('pyBigstick/nuclear_data.json')

# script comments
option_comment = '! create densities\n'
name_comment = '! output file name\n'
sps_comment = '! orbit information\n'
valence_comment = '! # of valence protons, neutrons\n'
jz2_comment = '! 2 * jz\n'
parity_comment = '! both parities wanted\n'
fragsize_comment = '! limit on fragment size for breaking Lanczos vectors\n'
hamil_comment = '! interaction file name (pkuo.int)\n'
scaling_comment = '! scaling\n'
end_comment = '! end of reading in Hamiltonian files\n'
diag_comment = '! diagonalization algorithm\n'
states_comment = '! number of state, max Lanczos iteration'


class Nucleus:
    def __init__(self, nucl_symbol: str, n_states=6, diag='ld', maxiter=400, fragsize=-1):
        symbol, mass_number = re.split('(\d+)', nucl_symbol)[:2]

        try:
            nucl = nuclear_data[nuclear_data['symbol']== symbol.capitalize()]
            proton = int(nucl['atomicNumber'])
            neutron = int(mass_number) - proton
            if neutron < 0:
                raise ValueError('Nucleus does not exist')
        except:
            raise ValueError('Nucleus does not exist')

        self.symbol = symbol
        self.nucl_symbol = nucl_symbol
        self.nucl = nucl
        self.A = int(mass_number)
        self.p = proton
        self.n = neutron

        # proton and neutron must be in the same orbit
        self.p_orbit = self.__get_orbit(proton)
        self.n_orbit = self.__get_orbit(neutron)
        if self.p_orbit != self.n_orbit:
            raise ValueError("Protons and neutrons are in different orbit. Handle them manually.")
        else:
            self.orbit = self.p_orbit

        self.jz = 0 if self.A %2 == 0 else 1
        self.p_valence = self.__get_valence()[0]
        self.n_valence = self.__get_valence()[1]

        self.int = self.__get_interaction()
        self.scaling = f'{1} {2*self.__get_core_orbit_capcity()+2} {self.A} {0.3}'

        self.n_states = n_states # number of states,
        self.diag = diag # diagonalize algo
        self.maxiter = maxiter #max iteration
        self.fragsize = fragsize # for parallel job only, -1 means no parallel

        # results
        self.states = [] # energy level states: staten, E, Ex, J, T
        self.densities = [] # density matrices: statei, statej, orba, orbb, Jt, Tt, value


    # n is nucleon
    def __get_valence(self, hole=True):
        core_size = self.__get_core_orbit_capcity()
        size = df_orbits.loc[df_orbits.name == self.orbit, 'size'].values[0]
        p_valence = self.p - core_size
        n_valence = self.n - core_size

        # negatives mean holes
        if hole == True:
            if p_valence > size/2:
                p_valence = - (size - p_valence)
            if n_valence > size/2:
                n_valence = - (size - n_valence)

        return p_valence, n_valence

    def __get_orbit(self, n):
        if n <= 2:
            return 's'
        if 2 < n <= 8:
            return '0p'
        if 8 < n <= 20:
            return 'sd'
        if 20 < n <= 40:
            return 'pf'
        if 50 < n <= 82:
            return 'jj55'
        else:
            raise ValueError('Too many nucleons')

    def __get_core_orbit_capcity(self):
        return df_orbits.loc[df_orbits.name == self.orbit, 'core_size'].values[0]

    def __get_interaction(self):
        if self.orbit == 's':
            raise ValueError('Interaction not found')
        return df_orbits.loc[df_orbits.name == self.orbit, 'int'].values[0]

    # generate the script
    def script_gen(self):
        option, end, parity = 'd', 'end', '0'
        output = f'{option:20} {option_comment}'
        output += f'{self.nucl_symbol:20} {name_comment}'
        output += f'{self.orbit:20} {sps_comment}'
        output += f'{str(self.p_valence)} {str(self.n_valence):18} {valence_comment}'
        output += f'{str(self.jz):20} {jz2_comment}'
        output += f'{parity:20} {parity_comment}'
        if self.fragsize != -1:
            output += f'{str(self.fragsize):20} {fragsize_comment}'
        if self.int == 'jj55pna':
            output += 'upn\n'
        output += f'{self.int:20} {hamil_comment}'
        output += f'{self.scaling:20} {scaling_comment}'
        output += f'{end:20} {end_comment}'
        output += f'{self.diag:20} {diag_comment}'
        output += f'{str(self.n_states)} {str(self.maxiter):18} {states_comment}'
        return output

    # save the script
    def script_save(self):
        filepath = f'examples/{self.nucl_symbol}/create_wfn.in'

        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            f.write(self.script_gen())

    # prepare necessary input files (.sps, .int, etc)
    def prepare(self):
        commands = f'cp sps/{self.orbit}.sps examples/{self.nucl_symbol}/;'
        commands += f'cp ints/{self.int}.int examples/{self.nucl_symbol}/;'
        os.system(commands)

    def run(self, bs=None, quiet=True):
        if bs == None:
            raise ValueError('Must specify BIGSTICK path by bs=/bs/path')
        commands = f'cd examples/{self.nucl_symbol};'
        commands += f'{bs} < create_wfn.in &&'
        if quiet == True:
            commands += f'rm *.bigstick fort.* {self.nucl_symbol}.lcoef;\n'
        os.system(commands)


    def __get_levels(self):
        filepath = f'examples/{self.nucl_symbol}/{self.nucl_symbol}.res'
        with open(filepath, "r") as f:
            lines = f.readlines()
            split_lines = [line.split() for line in lines]
            for line_n ,line in enumerate(split_lines):
                if 'State' in line and 'Ex' in line:
                    headings_n = line_n
                    break

            # the states
            for i in range(headings_n+1, headings_n+1+self.n_states):
                state = np.array(split_lines[i])
                state_n, energy, energy_x, J0, T0 = state
                self.states.append([int(state_n), float(energy), float(energy_x), float(J0), float(T0)])

        return self.states

    # return ij transition in bigstick format
    def __get_transition_ij(self, statej=2, statei=1):
        filepath = f'examples/{self.nucl_symbol}/{self.nucl_symbol}.dres'
        matrix = ''
        starting_line, ending_line = 0, 0
        with open(filepath) as f:
            unsplit_lines = f.readlines()
            lines = []
            for line in unsplit_lines:
                line = line.split()
                if line and '++++' in line[0]:
                    line = []
                lines.append(line)

            for i in range(len(lines)):
                if 'Initial' in lines[i] and 'state' in lines[i]:
                    if statei == int(lines[i][3]) and statej == int(lines[i+1][3]):
                        starting_line = i

            while True:
                line_pivot = starting_line + ending_line
                if lines[line_pivot] == []:
                    break
                ending_line += 1
            ending_line = starting_line + ending_line

            for l in range(starting_line, ending_line):
                matrix += unsplit_lines[l]
        return matrix

    def __get_matrices_ij(self, statej=2, statei=1):
        matrix = self.__get_transition_ij(statej, statei)
        lines = matrix.split('\n')
        lines = [line.split() for line in lines]

        for line in lines:
            if 'Jt' in line:
                Jt = int(line[2][:-1])
            if len(line) == 4:
                orba, orbb, me0, me1 = line
                self.densities.append([statei, statej, int(orba), int(orbb), Jt, 0, float(me0)])
                self.densities.append([statei, statej, int(orba), int(orbb), Jt, 1, float(me1)])

        return self.densities

    def save_results(self):
        # initilize
        self.states = []
        self.densities = []

        self.__get_levels()
        for i in range(self.n_states):
            for j in range(self.n_states):
                self.__get_matrices_ij(j, i)

        self.states = pd.DataFrame(self.states, columns=['state', 'E', 'Ex', 'J', 'T'])
        self.densities = pd.DataFrame(self.densities, columns=['statei', 'statej', 'orba', 'orbb', 'Jt', 'Tt', 'value'])


if __name__ == "__main__":
    nu = Nucleus('fe57')
    print('jz:', nu.jz)
