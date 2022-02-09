import re
import pandas as pd

nuclear_data = pd.read_json('./nuclear_data.json')

orbits = {
    # s
    "0s1/2": 2,

    # p
    "0p3/2": 4,
    "0p1/2": 2,

    # sd
    "0d5/2": 6,
    "1s1/2": 2,
    "0d3/2": 4,

    # pf
    "0f7/2": 8,
    "1p3/2": 4,
    "0f5/2": 6,
    "1p1/2": 2,

    "0g9/2": 10,

    # jj55
    "0g7/2": 8,
    "1d5/2": 6,
    "0h11/2": 12,
    "1d3/2": 4,
    "2s1/2": 2
}

# Input nucleon number
# Output the valence orbits
def orbit(n):
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


# Input element with mass number, eg Fe57, Ar40
# Output proton and neutron numbers
def nucleon(element: str):
    symbol, mass_number = re.split('(\d+)', element)[:2]
    nucl = nuclear_data[nuclear_data['symbol']== symbol.capitalize()]
    proton = int(nucl['atomicNumber'])
    neutron = int(mass_number) - proton
    return proton, neutron


if __name__ == "__main__":
    print(nucleon('fe57')[0])