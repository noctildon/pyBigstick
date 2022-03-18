import re
import pandas as pd
from orbits import df_orbits

nuclear_data = pd.read_json('./nuclear_data.json')

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
    def __init__(self, nucl_symbol: str, states=6, diag='ld', maxiter=400, fragsize=-1):
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
        self.A = mass_number
        self.p = proton
        self.n = neutron

        # proton and neutron must be in the same orbit
        self.p_orbit = self.get_orbit(proton)
        self.n_orbit = self.get_orbit(neutron)
        if self.p_orbit != self.n_orbit:
            raise ValueError("Protons and neutrons are in different orbit. Handle them manually.")
        else:
            self.orbit = self.p_orbit

        self.jz = self.get_jz()
        self.p_valence = self.get_valence()[0]
        self.n_valence = self.get_valence()[1]

        self.int = self.get_interaction()
        self.scaling = f'{1} {2*self.get_core_orbit_capcity()+2} {self.A} {0.3}'

        # number of states, diagonalize algo, max iteration
        self.states = states
        self.diag = diag
        self.maxiter = maxiter

        # for parallel job only, -1 means no parallel
        self.fragsize = fragsize

    # n is nucleon number
    def get_valence(self, hole=True):
        core_size = self.get_core_orbit_capcity()
        size = df_orbits.loc[df_orbits.name == self.orbit, 'size'].values[0]
        p_valence = self.p - core_size
        n_valence = self.n - core_size

        # negative valence means hole
        if hole == True:
            if p_valence > size/2:
                p_valence = - (size - p_valence)
            if n_valence > size/2:
                n_valence = - (size - n_valence)

        return p_valence, n_valence

    def get_orbit(self, n):
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

    def get_jz(self):
        if (self.p + self.n)%2 == 0:
            return 0
        else:
            return 1

    def get_core_orbit_capcity(self):
        return df_orbits.loc[df_orbits.name == self.orbit, 'core_size'].values[0]

    def get_interaction(self):
        if self.orbit == 's':
            raise ValueError('Interaction not found')
        return df_orbits.loc[df_orbits.name == self.orbit, 'int'].values[0]

    def script(self):
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
        output += f'{str(self.states)} {str(self.maxiter):18} {states_comment}'
        return output


if __name__ == "__main__":
    nu = Nucleus('fe57')
    print('p_valence:', nu.p_valence)
    print('n_valence:', nu.n_valence)
    print('jz:', nu.jz)
    print('scaling:',nu.scaling)
    print('**********')
