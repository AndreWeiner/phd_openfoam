Info << "Parallel consecutive reaction:" << endl;
Info << "A + B -(kAB)> P" << endl;
Info << "P + B -(kPB)> S" << endl;
Info << "A - Transfer species" << endl;
Info << "B - Bulk species" << endl;
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

Info<< "Reading scalar field B\n" << endl;

volScalarField B
(
    IOobject
    (
        "B",
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

Info<< "Reading scalar field S\n" << endl;

volScalarField S
(
    IOobject
    (
        "S",
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

Info<< "Reading diffusivity DB\n" << endl;

dimensionedScalar DB
(
    transportProperties.lookup("DB")
);

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    transportProperties.lookup("DP")
);

Info<< "Reading diffusivity DS\n" << endl;

dimensionedScalar DS
(
    transportProperties.lookup("DS")
);

Info<< "Reading reaction rate kAB\n" << endl;

dimensionedScalar kAB
(
    transportProperties.lookup("kAB")
);

Info<< "Reading reaction rate kAP\n" << endl;

dimensionedScalar kAP
(
    transportProperties.lookup("kAP")
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

const word B_model_name
(
    transportProperties.lookup("B_model")
);

const word S_model_name
(
    transportProperties.lookup("S_model")
);

const word rAB_model_name
(
    transportProperties.lookup("rAB_model")
);

const word rAP_model_name
(
    transportProperties.lookup("rAP_model")
);

torch::jit::script::Module A_model = torch::jit::load(A_model_name);
torch::jit::script::Module P_model = torch::jit::load(P_model_name);
torch::jit::script::Module B_model = torch::jit::load(B_model_name);
torch::jit::script::Module S_model = torch::jit::load(S_model_name);
torch::jit::script::Module rAB_model = torch::jit::load(rAB_model_name);
torch::jit::script::Module rAP_model = torch::jit::load(rAP_model_name);

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

surfaceScalarField B_Diff
(
    IOobject
    (
        "B_Diff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    DB
);

surfaceScalarField phiB
(
    IOobject
    (
        "phiB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

surfaceScalarField S_Diff
(
    IOobject
    (
        "S_Diff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    DS
);

surfaceScalarField phiS
(
    IOobject
    (
        "phiS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

torch::Tensor fTensorAB = torch::ones({patchFaceIDs.size(), 3}, torch::kFloat64);
torch::Tensor fTensorAP = torch::ones({patchFaceIDs.size(), 3}, torch::kFloat64);
torch::Tensor fTensorABP = torch::ones({patchFaceIDs.size(), 4}, torch::kFloat64);
torch::Tensor fTensorABPS = torch::ones({patchFaceIDs.size(), 5}, torch::kFloat64);
