# Create input for Bigstick
from helper import orbit, nucleon


option = 'd'
option_comment = '! create densities\n'

name = 'O15'
name_comment = '! output file name\n'



sps = '0p'
sps_comment = '! orbit information\n'
valence = '6 -1'
valence_comment = '! # of valence protons, neutrons\n'

p, n = nucleon(name)
sps_p = orbit(p)
sps_n = orbit(n)
if sps_p == sps_n:
    sps = sps_p
else:
    raise Exception("Protons and neutrons are in different orbit. Handle them manually.")


jz2 = '1'
jz2_comment = '! 2 * jz(spin of nucleus)\n'

hamil = 'pkuo'
hamil_comment = '! interaction file name (pkuo.int)\n'

scaling = '1 6 15 0.3' # p.44 in document
scaling_comment = '! scaling\n'

end = 'end'
end_comment = '! end of reading in Hamiltonian files\n'

diag = 'ex'
diag_comment = '! diagonalization algorithm\n'

states = '16'
states_comment = '! number of state'



# output = f'{option:20} {option_comment}'
# output += f'{name:20} {name_comment}'
# output += f'{sps:20} {sps_comment}'
# output += f'{valence:20} {valence_comment}'
# output += f'{jz2:20} {jz2_comment}'
# output += f'{hamil:20} {hamil_comment}'
# output += f'{scaling:20} {scaling_comment}'
# output += f'{end:20} {end_comment}'
# output += f'{diag:20} {diag_comment}'
# output += f'{states:20} {states_comment}'
# print(output)


