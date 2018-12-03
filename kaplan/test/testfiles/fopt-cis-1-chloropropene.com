%mem=2GB
%chk=/work/jgarne01/assign2/cis-1-chloropropene/B3LYP/cis-1-chloropropene-B3LYP-Def2TZVP.chk
%nproc=8
#B3LYP/Def2TZVP opt

notes
cis-1-chloropropene

0 1
Cl  -1.6393    0.6138    0.0000
C    1.4577    0.7156    0.0000
C    0.7556   -0.6017    0.0000
C   -0.5740   -0.7278    0.0000
H    0.8089    1.5936    0.0000
H    2.1004    0.7829    0.8841
H    2.1003    0.7829   -0.8841
H    1.3696   -1.4976    0.0000
H   -1.0898   -1.6783   -0.0001

