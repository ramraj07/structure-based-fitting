# structure-based-fitting
Demonstration of structure-based superresolution data refinement method used in the Single Molecule Localization Challenge

This repository contains a Matlab Live Notebook, and HTML/PDF versions of the Live Notebook, 
which explain the method by which the [SMLM challenge](http://bigwww.epfl.ch/smlm/challenge2016/) was solved for 3D data for our lab submission.

The GIST is that I used k-means clustering and manual classification of clusters as belonging to different unique biological structures 
to identify which single molecule localizations belong to the same structure. Then it was easy to use location information 
from each single molecule to enhance the accuracy of the location estimate of neaighboring molecules. 

Given that the "axial" (Z) dimension estimation accuracy is always an order of magnitude less accurate than the lateral dimension 
accuracies, this method allowed the "transfer" of the accuracy gained in the lateral dimensions to the axial dimension. 
