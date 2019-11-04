Info << "Parallel consecutive reaction:" << endl;
Info << "A + B -(kAB)> P" << endl;
Info << "P + B -(kPB)> S" << endl;
Info << "A - Transfer species" << endl;
Info << "B - Bulk species" << endl;
Info << "Reading parallelConsecutiveReaction properties." << endl;

IOdictionary parallelConsecutiveReaction
(
    IOobject
    (
        "parallelConsecutiveReaction",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

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

Info<< "Reading scalar field Aphy\n" << endl;

volScalarField Aphy
(
    IOobject
    (
        "Aphy",
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
    parallelConsecutiveReaction.lookup("DA")
);

Info<< "Reading diffusivity DB\n" << endl;

dimensionedScalar DB
(
    parallelConsecutiveReaction.lookup("DB")
);

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    parallelConsecutiveReaction.lookup("DP")
);

Info<< "Reading diffusivity DS\n" << endl;

dimensionedScalar DS
(
    parallelConsecutiveReaction.lookup("DS")
);

Info<< "Reading reaction rate kAB\n" << endl;

dimensionedScalar kAB
(
    parallelConsecutiveReaction.lookup("kAB")
);

Info<< "Reading reaction rate kAP\n" << endl;

dimensionedScalar kAP
(
    parallelConsecutiveReaction.lookup("kAP")
);
