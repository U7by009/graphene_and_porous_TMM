# graphene_and_porous_TMM
Matlab functions for Transfer Matrix Method using graphene and a porous layer.

Matlab functions to visualise the reflectance in a multilayer dielectric structure, where the second layer consists of a stack of a number of graphene monolayers and the penultimate layer is a porous dielectric sheet. The porosity of the dielectric medium is evaluated using Bruggerman's formula. The input parameters are:

The number N of graphene monolayers (modelled as a dielectric).

The porosity P of the penultimate layer in contact with the sensing aquous material.

The operating wavelength in nanometres.

The chemical potential (Fermi energy) in meV.

The temperature of the system \textbf{T}, measured in K.

The relaxation frequency of the system \textbf{gamma} in rad-s$^{-1}$.

A vector of the refractive indices of the materials in each layer [indices(1) indices(2) ... indices(end)]. The index of the porous layer (end-1) must be written as the corresponding index for a layer without pores (P=0).

A vector with the thicknesses of each layer [layers(1) layers(2) ... layers(end)] specified in nanometres. The thicknesses of layers(1) and layers(end) must be entered as NaN (Not-a-number) since they are considered semi-infinite layers.

A vector with the angles of incidence at the first interface [theta(1) theta(2) ... theta(end)] specified in radians.

# Matlab functions
TMM.m: This code visualises the dependence of the reflectance for both TE and TMM polarised light as a function of the angle of incidence for a given porosity of the penultimate layer. For this case, graphene is considered as a dielectric medium (n+ikappa) of thickness N-0.34 nm, without a current density in the dielectric interlayer (sigma=0).

TMMcontour.m: This code visualises the dependence of the reflectance for both TE and TMM polarised light as a function of the angle of incidence and porosity. For this case, graphene is considered as a dielectric medium (n+ikappa) of thickness N-0.34 nm, without a current density in the dielectric interlayer (sigma=0).

TMM2.m: This code visualises the dependence of the reflectance, both for TE and TMM polarised light, as a function of the angle of incidence, for a given porosity of the penultimate layer. For this case, graphene is considered as a dielectric medium with no extinction coefficient (epsilon'=n^2, epsilon''=0) of thickness N-0.34 nm, with a current density in the dielectric interlayer (sigma~=0) given by the Kubo model for the inter-bands and intra-band phenomena in graphene.

TMM2contour.m: This code visualises the dependence of the reflectance for both TE and TMM polarised light as a function of the angle of incidence and porosity.  For this case, graphene is considered as a dielectric medium with no extinction coefficient (epsilon'=n^2, epsilon''=0) of thickness N-0.34 nm, with a current density in the dielectric interlayer (sigma~=0) given by the Kubo model for the inter-bands and intra-band phenomena in graphene.

TMM3.m: This code visualises the dependence of the reflectance for both TE and TMM polarised light as a function of the angle of incidence for a given porosity of the penultimate layer. For this case, graphene is considered jointly as a dielectric medium (epsilon'+epsilon''i) of thickness N-0.34 nm and with a current density in the dielectric inter-faces (sigma~=0) given by the Kubo model for the inter-bands and intra-band phenomena in graphene.

TMM3contour.m: This code visualises the dependence of the reflectance for both TE and TMM polarised light as a function of the angle of incidence and porosity. For this case, graphene is considered jointly as a dielectric medium (epsilon'+epsilon''i) of thickness N-0.34 nm and with a current density in the dielectric inter-faces (sigma~=0) given by the Kubo model for the inter-bands and intra-band phenomena in graphene.
