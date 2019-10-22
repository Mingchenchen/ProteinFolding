## Introduction

A protein can be simplified to a string of beads, which either correspond to a hydrophobic (H) or a polar (P) amino acid residue. Any pair of beads can interact by a [Lennard-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) V(r) given by:

<p align="center">V(r) = e[r<sup>-12</sup>-2r<sup>-6</sup>]</p>

where r is the distance between the two beads. Assuming both beads are not direct neighbours in the primary structure, we let:
* e = E<sup>HH</sup> = -2.3 if both beads are of type H.
* e = E<sup>HP</sup> = -1 if both beads are of different types.
* e = E<sup>PP</sup> = 0 if both beads are of type P.

Otherwise, if they are direct neighbours, we let e = E<sup>NN</sup> = -10 independent of type.

## Methods

To simulate protein folding in 2-D, we implemented a Metropolis algorithm, which accepts a trial move (small change in the position of a randomly chosen bead in the x-y plane) with a probability:

<p align="center">p('accept') = min{1, exp(-(E<sub>f</sub>-E<sub>i</sub>)/T)}</p>,

where E<sub>f</sub> is the new energy after the trial move and E<sub>i</sub> is the initial energy before the trial move. T is the temperature. To simplify our model, we initialise T as 1, and consider only the bead that moves when calculating the energy change during a trial move.

To test our code, we used the following 36-residue long polypeptide [(Li et al., 1996)][https://science.sciencemag.org/content/273/5275/666.long]:

## Credits

This repository was inspired by a systems biology coursework from the final year of my bachelor's studies in biotechnology, majoring in computational biology. During this course, I learnt programming from scratch and, as a rookie, I was not able to solve this problem at all back then. The lecturer of this course was [Prof. Robert Endres](https://rgendres3.wixsite.com/biologicalphysics).

The image of the test polypeptide was taken from a paper by Ned Wingreen's group. The full reference of this paper is:
* Li, H., Helling, R., Tang, C., & Wingreen, N. (1996). Emergence of Preferred Structures in a Simple Model of Protein Folding. Science, 273(5275), 666–669.
