## Introduction

A protein can be simplified to a string of beads, which either correspond to a hydrophobic (H) or a polar (P) amino acid residue. Any pair of beads can interact by a [Lennard-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) V(r) given by:

<p align="center">V(r) = e[r<sup>-12</sup>-2r<sup>-6</sup>]</p>

where r is the distance between the two beads. Assuming both beads are not direct neighbours in the primary structure, we let:
* e = E<sup>HH</sup> = -2.3 if both beads are of type H.
* e = E<sup>HP</sup> = -1 if both beads are of different types.
* e = E<sup>PP</sup> = 0 if both beads are of type P.

Otherwise, if they are direct neighbours, we let e = E<sup>NN</sup> = -10 independent of type.

By summing up the Lennard-Jones potential of all pairs of beads, we can calculate the difference in the energy before and after a given folding sequence. This allows us to show that protein folding is driven by a decrease in energy (a.k.a. a release of energy to the surrounding environment).

## Why is protein folding considered an unsolved problem in biology?



## Methods

The code for simulating simulate protein folding in 2-D is found in [`prot_folding_NSEW.m`](https://github.com/liweiyap/ProteinFolding/blob/master/prot_folding_NSEW.m). We implemented the Metropolis algorithm, which accepts a trial move (small change in the position of a randomly chosen bead in the x-y plane) with a probability:

<p align="center">p('accept') = min{1, exp(-(E<sub>f</sub>-E<sub>i</sub>)/T)}</p>

where E<sub>f</sub> is the new energy after the trial move and E<sub>i</sub> is the initial energy before the trial move. T is the temperature. To simplify our model, we initialise T as 1, and consider only the bead that moves when calculating the energy change during a trial move. The number of trial moves, i.e. the number of iterations of the Metropolis algorithm, is user-defined.

To test our code, we used the following 36-residue long polypeptide [(Li et al., 1996)](https://science.sciencemag.org/content/273/5275/666.long):

![alt text](https://github.com/liweiyap/ProteinFolding/blob/master/test_polypeptide.png)

A sample output of the code following 10<sup>7</sup> iterations of the Metropolis algorithm is shown [here](https://github.com/liweiyap/ProteinFolding/blob/master/sample_output_with_ten_million_steps.png). The computational time taken was 145 seconds. In a folded protein, one would expect the hydrophobic residues (in red) to be shielded from the aqueous environment by the polar residues (in blue). However, it is possible that we do not get such a perfect folding because our polypeptide in question is only 36 residues long. On the other hand, hydrophobic patches have been reported to exist on the external surface of proteins and might even have important biological functions.

## Credits

This repository was inspired by a systems biology coursework from the final year of my bachelor's studies in biotechnology, majoring in computational biology. During this course, I learnt programming from scratch and, as a rookie, I was not able to solve this problem at all back then. The lecturer of this course was [Prof. Robert Endres](https://rgendres3.wixsite.com/biologicalphysics).

The image of the test polypeptide was taken from a paper by Ned Wingreen's group. The full reference of this paper is:
* Li, H., Helling, R., Tang, C., & Wingreen, N. (1996). Emergence of Preferred Structures in a Simple Model of Protein Folding. Science, 273(5275), 666–669.
