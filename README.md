# Heteroclinics_PtoP_1D
Matlab codes to reproduce results in "Heteroclinic connections between finite-amplitude periodic orbits emerging from a codimension two singularity" by T.J. Bridges, D.J.B. Lloyd, D.J. Ratliff, and P. Sprenger. ArXiv link to preprint: https://arxiv.org/abs/2602.05680

Numerical continuation codes use the following repository
- Daniele Avitabile’s continuation-spatially-extended-systems-tutorial: DOI 10.5281/zenodo.3821163 - https://zenodo.org/records/3821169

1. SHE357_foliations directory: Computes the unstable and stable manifolds of the periodic orbits in the Swift-Hohenberg equation reproducing Figure 5 in the paper
- To run: run script > run_compute_foliations.m

2. SHE357_codim2_pts directory: Computes the Hamiltonian codim 2 points of the periodic orbits in the Swift-Hohenberg equation reproducing Figure 12 in the paper
- To run you need to run the following scripts sequentially:
- run_SHE357_periodic.m
- run_SHE357_Hk_0.m
- run_SH357_codim2_point.m
There is a detailed readme file in the directory

3. Boussinesq Fronts directory: Computes the Periodic-to-Periodic heteroclinics in the AC-Boussinesq equations reproducing Figure 16 in the paper. The codes are set up for general two-coupled second order reaction-diffusion systems and can be adapted to other equations e.g. SHE357.
- To run: run script > run_get_PtoP_front_ac_Boussinesq.m
