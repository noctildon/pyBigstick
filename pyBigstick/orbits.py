from numpy import NaN
import pandas as pd

orbits = {
    'name': ['s', '0p', 'sd', 'pf', 'pfg', 'jj55'],
    'size': [2, 6, 12, 20, 22, 32],
    'core_size': [0, 2, 8, 20, 28, 50],
    'int': [NaN, 'pkuo', 'usdb', 'GXPF1', 'JUN45', 'jj55pna']
}

df_orbits = pd.DataFrame(data=orbits)


if __name__ == "__main__":
    print(df_orbits.loc[df_orbits.name == 'sd', 'core_size'].values[0])