# sgsReactingPimpleFoam
## Overview
sgsReactingPimpleFoam is a pimpleFoam-based solver for the incompressible
Navier-Stokes equations with additional transport of chemical species. The species
transport is currently hardcoded (instead of using the OpenFOAM chemistry infrastructure)
to facilitate the development and testing of subgrid-scale (SGS) models for convection-dominated
boundary layers. The solver at the moment of writing only supports first-order decay reactions
and uses PyTorch to include machine learning models.

## General workflow

The SGS-related code is included as follows:
- initialize ML models (load models from *mlModels* folder)
- initialize SGS model (find interface cells and neighbors, compute important geometric quantities)
- correct convective and diffusive fluxes in interface cells via *D* and *phi* fields
- solve species transport for specified reaction
- write out relevant interfacial data (Sh number, correction factors, ...)

## Relevant files

- **sgsReactingPimpleFoam.C**: main solver file with *main* routine
- **initMLModel.H**: create and load ML models, run test forward pass through model, load ML model properties (min and max values for rescaling)
- **initSGSModel.H**: find interface cells, neighbors, neighboring faces
- **createFields.H**: read species properties, read SGS properties, create switches, allocate fields for post-processing output
- **createSpeciesFields.H**: read reaction and transport properties, create species fields
- **sgsFluxCorrection.H**: create fields for corrected transport properties (fluxes), compute and apply fluxes
- **solveSpeciesTransport.H**: solve reaction-dependent species transport equations, bound fields
- **printInterfaceData.H**: write out interfacial data regarding reactive mass transfer (local and global Sherwood number, polar angle, SGS correction factors, ...)
