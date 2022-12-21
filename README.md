# EKI-for-MRE
Ensemble Kalman Inversion (EKI) code for Magnetic Resonance Elastography (MRE)
Here we provide some of the code that we employed to produce the 
2D synthetic examples of the paper:
Marco Iglesias et al 2022 Phys. Med. Biol. 67 235003
that can be found:
https://iopscience.iop.org/article/10.1088/1361-6560/ac9fa1

This is research code distributed for the sake of clarity, transparency and 
data reproducibility of the above paper. The code is not optimised or 
customised for practical MRE applications. 

The main script is Driver_MRE.m. Running this script will:
(A) produce noisy displacements (Figure 4 in paper) from the true stiffness (saved in Truth.mat)
(B) Define prior and produce prior ensemble.
(C) Run EKI and plot the truth, posterior/prior mean/variance for Storage/Loss modulus (Figure 3 and Figure 6 in paper).

Most parameters are defined in Driver_MRE.m. but some parameters (assumed known) for the visco-elastic model are hard coded in set_model.m

Note that the main results from Section 4.2 are obtained using ensemble size N_En=10^4 which required us to use multiple cores from a node in the HPC at Nottingham. However, the code can still provide reasonable outcomes using a moderate size (i.e. N_En=500) more manageable for standard multi-core desktop/laptop.

The files tight_subplot.m and imagescwithnan.m were found online and used for visualisation. 
