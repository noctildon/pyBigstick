# Initialization

```py
nu = Nucleus(nucl_symbol='f19', n_state=6, diag='ld', maxiter=400, fragsize=-1)
```


# Method

```py
nu.script_save()  # generate BIGSTICK script and save in examples/
nu.prepare()      # copy the inputs
```

# Run

```py
bs = os.getcwd()+'/src/bigstick.x'    # specify the path
nu.run(bs=bs)                         # run the script
```