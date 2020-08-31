# cgtm
coarse grained transition matrix

Included: Sample pdb, processHop.out, and rts_fs.txt output. Not included: Sample trajectory. 

First: Voronoi 
You need:
1.	Voronoi package
2.	A pdb of your system containing only the membrane and receptor or receptor complex, stripped of nonpolar hydrogens. Segname of membrane should be “MEMB” and residue numbering should start at 1 and increase sequentially. “mypdb.pdb”
3.	A binary trajectory containing the same atoms as the pdb. “mytrajectory.dcd”

Steps:
1.	Compile Voronoi code. In hopBorders.C, edit the lipidType list to include the resnames of your lipid types, and return distinct numbers (0, 1, 2, etc.) for all residues you wish to have distinct grouping. For example, we assigned 0 to DPPC, 1 to Cholesterol, and 2 to DOPC. The numbers correspond to color codes.
2.	hopBorders mypdb.pdb mytrajectory.dcd > hopBorders.out
This produces “borders.bin” needed for step 3.
3.	processHop borders.bin 2 > processHop.out
The output file processHop.out contains a list of shell number and residue number for each lipid, listed frame-by-frame. Each frame is listed in a single line of text, followed by two summary statements which sum the lipids of each type in the first and second shell. 
Here is an example of the first 3 lines of the file:
7  1 4  2 2  3 4  4 6  5 3  6 6  7 6  8 1  9 4  10 4  11 4  12 5  13 4  14 2  15 2  16 1  17 3  18 6  19 3  20 4  21 5  22 4  23 5  24 5  25 2  26 5  27 8  28 6  29 6  30 5  31 8  32 2  33 …
    First shell: 15 15 9
Second shell: 30 17 7
In the first line, “7	1” indicates lipid with residue number 1 is in the 7th shell. “7	1  4	2  2	3” indicates lipids with residue numbers 1, 2, and 3 are in the 7th, 4th, and 2nd shells, respectively. “First shell: 15 15 9” counts the number of DPPC, Cholesterol, and DOPC molecules in the first shell. 
4.	python textpro_voronoi.py 
This produces “rts_fs.txt”, which contains a list of first shell and second shell tallies, with each line corresponding to a frame. The above example processHop.out would yield “15 15 9 30 17 7” in line one for tallies of DPPC, DOPC, and Cholesterol in the 1st and 2nd shells in the first frame. 

Second: CGTM
1.	Run the two cgtm codes, “cgtm1.m” and “cgtm2.m” sequentially in Matlab. This will compute the raw occupation probabilities and equilibrium distribution associated with the transitions from the input trajectory.
2.	To create a CGTM from multiple short trajectories, you will need to write a workflow that performs cgtm1.m on all the short trajectories and allows cgtm2.m to populate the transition matrix from your multiple input trajectories. This step in cgtm2.m should be operated on all trajectories, with the full state space accounted for:
    for i = j:n-1        
        oldstate = states_in_ts(i);
        newstate = difference(i) + states_in_ts(i);
        trans_prob(oldstate,newstate) = trans_prob(oldstate,newstate) + 1;
    end
