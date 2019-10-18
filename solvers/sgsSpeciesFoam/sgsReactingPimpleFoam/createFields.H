Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"

#include "createSpeciesFields.H"

IOdictionary speciesProperties
(
    IOobject
    (
        "speciesProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar refLength
(
    speciesProperties.lookup("referenceLength")
);

dimensionedScalar riseVelocity
(
    speciesProperties.lookup("riseVelocity")
);

scalar pi = constant::mathematical::pi;

Switch applySGS(speciesProperties.lookup("applySGSmodel"));

Switch diffSigmaCorr(speciesProperties.lookup("sgsDiffusionAtSigmaCorr"));
Switch diffFacesCorr(speciesProperties.lookup("sgsDiffusionFacesCorr"));
Switch advCorr(speciesProperties.lookup("sgsAdvectionCorr"));
Switch blendAdvection(speciesProperties.lookup("blendAdvection"));
Switch blendDiffBulk(speciesProperties.lookup("blendDiffBulk"));
Switch blendDiffSig(speciesProperties.lookup("blendDiffSig"));

dimensionedScalar alphaAdv
(
    speciesProperties.lookup("alphaAdv")
);
dimensionedScalar alphaDiff
(
    speciesProperties.lookup("alphaDiff")
);
dimensionedScalar alphaDiffSigma
(
    speciesProperties.lookup("alphaDiffSigma")
);

// Detect the free surface patch
label surfaceID(-1);
forAll (mesh.boundary(), patchI)
{
    if
    (
        (mesh.boundary()[patchI].name() == "freeSurface")
     || (mesh.boundary()[patchI].name() == "interface")
    )
    {
        surfaceID = patchI;
    }
}

if(surfaceID == -1)
{
    FatalErrorIn("reactingPimpleFoam")
        << "Free surface / interface patch not defined.  Please make sure that "
            << " the free surface patch is named as freeSurface or interface"
            << abort(FatalError);
}

// number of interface cells
int nInterfaceCells = fields[0].boundaryField()[surfaceID].size();

// output of model predictions
scalarField cBulkFinal(nInterfaceCells, 0.0);
scalarField dcBulkFinal(nInterfaceCells, 0.0);
scalarField dcSigFinal(nInterfaceCells, 0.0);
scalarField corrAdv(nInterfaceCells, 1.0);
scalarField corrFace(nInterfaceCells, 1.0);
scalarField corrSig(nInterfaceCells, 1.0);
