# import the library
from pyBigstick.nucleus import Nucleus

# This is BIGSTICK path
# Download BIGSTICK source code and compile it
bs = '/home/wei-chih/Desktop/Research/BigstickPublick/v7.9.12/bigstick.x'


# create nucleus, here f is fluorine (case insensitive), 19 is the mass number
nu = Nucleus('f19')

nu.script_save() # generate BIGSTICK script and save in examples/
nu.prepare() # copy the inputs
nu.run(bs=bs) # run the script