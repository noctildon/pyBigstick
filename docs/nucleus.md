# Initialization

```py
nu = Nucleus(nucl_symbol='f19', n_state=6, diag='ld', maxiter=400, fragsize=-1)
```

- nucl_symbol : nuclear symbol (or chemical symbol) that represents a unique nucleus type, followed by the mass number.
- n_state : the number of states that BIGSTICK will calculate
- diag : diagonalization method. 'ld' should be sufficient
- maxiter : the maximum iteration for each state
- fragsize : the fragment size in parallel job. set this to -1 for serial job


# Preparation

```py
nu.script_gen()   # return BIGSTICK input script in string
nu.script_save()  # save the input script to examples/f19
nu.prepare()      # copy other necessary input files
```

# Run

```py
# run the script. bs specifies the path, quiet turns off most of the message
nu.run(bs='/path/to/bigstick.x', quiet=True)
```
Then BIGSTICK will generate a few files, such as f19.dres, f19.res, f19.wfn. In which f19.dres has the most useful contents. Ignore all others.


# Complete
```py
nu.save_results()     # load results (f19.dres) to ram
print(nu.states)      # print the energy levels
print(nu.densities)   # print the density matrices
```


# Misc
```py
nu.check()   # check if the calculation has been done before
nu.clean()   # clean the previous output files
```