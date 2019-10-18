Info << "Single reaction:" << endl;
Info << "A + B -(kAB)> P" << endl;
Info << "A - Transfer species" << endl;
Info << "B - Bulk species" << endl;
Info << "Reading singleReaction properties." << endl;

IOdictionary singleReaction
(
    IOobject
    (
        "singleReaction",
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

Info<< "Reading diffusivity DA\n" << endl;

dimensionedScalar DA
(
    singleReaction.lookup("DA")
);

Info<< "Reading diffusivity DB\n" << endl;

dimensionedScalar DB
(
    singleReaction.lookup("DB")
);

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    singleReaction.lookup("DP")
);

Info<< "Reading reaction rate kAB\n" << endl;

dimensionedScalar kAB
(
    singleReaction.lookup("kAB")
);
