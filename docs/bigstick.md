### This is a introductory manual for BIGSTICK


# Inputs

## script
The script below provides options line by line to BIGSTICK
```
! exclamation mark followed by comments
! create_wfn.in
d             ! create density matrices (the .dres files)
f19           ! output name
sd            ! orbit information (sd.sps)
1 2           ! # of valence protons, neutrons (compare to O16)
1             ! 2 * jz(spin of nucleus), 0 for even mass number, 1 for odd mass number
0             ! both parities wanted
usdb          ! interaction file name (usdb.int)
1 18. 19. 0.3 ! scaling, 1 for single particle energy, 18 for core mass number+2, 19 for target mass (p.44 in document)
end           ! end of reading in Hamiltonian files
ld            ! Lanczos diagonalization algorithm
6 400         ! 6 states kept (max # iterations for lanczos)
```


```sh
# Dump the script to BIGSTICK
$ bigstick.x < create_wfn.in
```


## Nuclear files

### Particle space
Ending in .sps, stored in sps/, used to describe the orbits that nucleons are allowed to stay.


### Interaction
Ending in .int, stored in int/, used to describe the interaction of the nucleons in echo of all possible orbits.


# Outputs

## dres
In .dres files you will see two patterns, states and matrices

### States
```
  State      E        Ex         J       T
    1    -23.86096   0.00000     0.500   0.500
    2    -23.78367   0.07729     2.500   0.500
    3    -22.09059   1.77037     1.500   0.500
    4    -21.26237   2.59859     4.500   0.500
    5    -19.25724   4.60372     6.500   0.500
    6    -18.99335   4.86761     3.500   0.500
```

- State : the number of the state, 1 is the ground state
- E : the energy of the state in MeV
- Ex : the excitation energy in MeV
- J : the spin of the state
- T : the isospin of the state


### Matrices
```
 Initial state #    1 E =  -23.86096 2xJ, 2xT =    1   1
 Final state   #    1 E =  -23.86096 2xJ, 2xT =    1   1
 Jt =   0, Tt = 0        1
    1    1   0.19126   0.02824
    2    2   0.89604   0.34888
    3    3   1.17753   0.35578
 Jt =   1, Tt = 0        1
    1    1  -0.05634   0.01600
    1    2   0.13395  -0.10506
    1    3  -0.01476  -0.01563
    2    1  -0.13395   0.10506
    2    2   0.12703  -0.24666
    3    1   0.01476   0.01563
    3    3   0.42662  -0.39164
```
The 2D matrices writes the transition amplitude of 2 orbits for a given Jt, Tt, initial state and final state. For example, the second last line
```
3    1   0.01476   0.01563
```
This tells that, for Tt=0, the amplitude from orbit 3 to orbit 1 is 0.01476, and for Tt=1, the amplitude is 0.01563.

Jt and Tt are the spin and isospin of the transition operator, respectively.

One can use these matrices to compute the amplitude of any operator.