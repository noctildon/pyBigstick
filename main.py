import os
from pyBigstick.nucleus import Nucleus

# This is BIGSTICK path
# Copy BIGSTICK source code and compile it
bs = os.getcwd() +'/src/bigstick.x'


# Create nucleus, here f is fluorine (case insensitive), 19 is the mass number
nu = Nucleus('f19')


# Running BIGSTICK
nu.script_save() # generate BIGSTICK script and save in examples/
nu.prepare() # copy the inputs
nu.run(bs=bs) # run the script


# Save results to ram
nu.save_results()
print(nu.states) # print the energy levels
print('**************')
print(nu.densities) # print the density matrices
