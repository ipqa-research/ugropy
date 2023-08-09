# ugropy

TODO:

- Test all the combinations of the pyridine groups
- Add R, Q and main group to the unifac_smarts.csv
- More complex algorithm to handle the composed structure problem. (test set as "composed") 

Remember:
CHO it's for aldehyde
CH-O ether

Diferencias entre UNIFAC y PSRK:

[CH2]=[CH2] not in UNIFAC but yes in PSRK

https://opsin.ch.cam.ac.uk/

https://pubchem.ncbi.nlm.nih.gov//edit3/index.html

Ambiguos
test_13_ch2o.py  
("Benzyl 2-hydroxyethyl carbonate", {"CH2": 1, "OH": 1, "ACCH2": 1, "ACH": 5, "COO": 1, "CH2O": 1})

These problematic structures must be handled differently
[cH0][CH2]O[CX4H0]|"{""ACCH2"": -1, ""AC"": 1, ""CH2"": 1}"
[cH0][CH]O[CX4H0]|"{""ACCH"": -1, ""AC"": 1, ""CH"": 1}"

The difference between a normal problematic structure is that the problem of 
a normal problematic structure relies on the impossibility to construct short 
and simple SMARTS that represent the functional groups and doesn't create an 
impossible-to-handle false positive by the general logic of the algorithm.

On the other hand, the above structures "can be handled" by the algorithm, the
problem in that specific case is that the ether group CH0-O doesn't exist in 
the classical UNIFAC model. By that, the algorithm has to break its own rule
of using the most specific functional groups always. The most specific group
ACCH2 cannot be used and has to be separated into the less specific groups: AC 
and CH2. Since the algorithm will also detect the other specific group CH2O
the replacement of the ACCH2 by an AC CH2 (decomposing the functional group) 
will do the thing.

The specific groups like ACCH2 will be called composed groups. Understanding
a composed group as a functional group that can be represented by two less 
specific groups, e.g:

- ACCH2 -> AC + CH2
- ACOH -> AC + OH
- CH3COO -> CH3 + COO
