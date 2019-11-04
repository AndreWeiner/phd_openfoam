Info << "Decay reaction:" << endl;
Info << "A -(kA)> P" << endl;
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

Info<< "Reading scalar field P\n" << endl;

volScalarField P
(
    IOobject
    (
        "P",
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

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    transportProperties.lookup("DP")
);

Info<< "Reading reaction rate kA\n" << endl;

dimensionedScalar kA
(
    transportProperties.lookup("kA")
);

Info<< "Reading subgrid-scale model(s)\n" << endl;

const word A_model_name
(
    transportProperties.lookup("A_model")
);

const word P_model_name
(
    transportProperties.lookup("P_model")
);

torch::jit::script::Module A_model = torch::jit::load(A_model_name);
torch::jit::script::Module P_model = torch::jit::load(P_model_name);

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

surfaceScalarField P_Diff
(
    IOobject
    (
        "P_Diff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    DP
);

surfaceScalarField phiP
(
    IOobject
    (
        "phiP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

torch::Tensor fTensorA = torch::ones({patchFaceIDs.size(), 2}, torch::kFloat64);
torch::Tensor fTensorAP = torch::ones({patchFaceIDs.size(), 3}, torch::kFloat64);
