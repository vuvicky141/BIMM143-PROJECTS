new12
================
Vicky Vu
2/13/2020

## PDB Statistics

Here we inspect the types of structures in themain database for 3D
biomolecular data - the PDB.

Q1: Determine the percentage of structures solved by X-Ray and Electron
Microscopy. Also can you determine what proportion of structures are
protein?

``` r
#read in the file 
stats <- read.csv("Data Export Summary.csv", row.names = 1)
stats
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 133756          2086               6884     8 142734
    ## NMR                    11308          1317                265     8  12898
    ## Electron Microscopy     3241            35               1095     0   4371
    ## Other                    284             4                  6    13    307
    ## Multi Method             146             5                  2     1    154

``` r
#calculate percentage values asked for  
sum  <- sum (stats [1,5],stats[3,5])
total <- sum ( stats [,5])

sum/total *100
```

    ## [1] 91.67477

``` r
# percentage of each 
stats$Total/sum(stats$Total)*100
```

    ## [1] 88.95079270  8.03793997  2.72397547  0.19132017  0.09597168

``` r
# what porportion of structures are protein 
sum(stats$Proteins)/sum(stats$Total)*100
```

    ## [1] 92.69057

Q2: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?
1289 HIV-1 protease structures found in te current PDB.

# part 2

Read a single PDB structure into R

``` r
library(bio3d)

pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

# Select protein only

Heere I will use the ‘atom.select()’ function to do this…

``` r
lig <- atom.select(pdb, "ligand", value = TRUE)

prot <- atom.select(pdb, "protein", value = TRUE)
```

…and write out these new PDB objects with the "write.pdb()’ function.

``` r
write.pdb(lig, file = "1hsg_ligand.pdb")
write.pdb( prot, file = "1hsg_protein.pdb")
```

``` r
attributes(pdb)
```

    ## $names
    ## [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  
    ## 
    ## $class
    ## [1] "pdb" "sse"
