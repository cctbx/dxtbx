dxtbx Experiment: add the concept of an index property, which can be used as the
index in an experiment list, such that if the list is split the correct 
reflections can be matched up with that experiment. 

N.B. implication of this is to replace this pattern:

```
for j, expt in enumerate(experiments):
    refl = refl.select(refl["id"] == j)
```

with 

```
for expt in experiments:
    refl = refl.select(refl["id"] == expt.index)
```

N.B. also that this property is ephermeral i.e. is only valid for the run-time 
of the program and is not serialised to disk (where experiment identifiers 
should be used)

