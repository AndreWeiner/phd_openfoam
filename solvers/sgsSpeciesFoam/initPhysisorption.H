Info << "Physisorption:" << endl;
Info << "A - Transfer species" << endl;
Info << "Reading transport and reaction properties." << endl;

Info<< "Reading scalar field A\n" << endl;

volScalarField A
(
    IOobject
    (
        "A",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading diffusivity DA\n" << endl;

dimensionedScalar DA
(
    transportProperties.lookup("DA")
);

Info<< "Reading subgrid-scale model(s)\n" << endl;

const word A_model_name
(
    transportProperties.lookup("A_model")
);

torch::jit::script::Module A_model = torch::jit::load(A_model_name);

surfaceScalarField A_Diff
(
    IOobject
    (
        "A_Diff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    DA
);

surfaceScalarField phiA
(
    IOobject
    (
        "phiA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

torch::Tensor fTensorA = torch::ones({patchFaceIDs.size(), 2}, torch::kFloat64);
