Info << "Parallel competitive reaction:" << endl;
Info << "A + B -(kAB)> P" << endl;
Info << "A + C -(kAC)> S" << endl;
Info << "A - Transfer species" << endl;
Info << "B/C - Bulk species" << endl;
Info << "Reading parallelCompetitiveReaction properties." << endl;

IOdictionary parallelCompetitiveReaction
(
    IOobject
    (
        "parallelCompetitiveReaction",
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

Info<< "Reading scalar field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
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
    parallelCompetitiveReaction.lookup("DA")
);

Info<< "Reading diffusivity DB\n" << endl;

dimensionedScalar DB
(
    parallelCompetitiveReaction.lookup("DB")
);

Info<< "Reading diffusivity DC\n" << endl;

dimensionedScalar DC
(
    parallelCompetitiveReaction.lookup("DC")
);

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    parallelCompetitiveReaction.lookup("DP")
);

Info<< "Reading diffusivity DS\n" << endl;

dimensionedScalar DS
(
    parallelCompetitiveReaction.lookup("DS")
);

Info<< "Reading reaction rate kAB\n" << endl;

dimensionedScalar kAB
(
    parallelCompetitiveReaction.lookup("kAB")
);

Info<< "Reading reaction rate kAC\n" << endl;

dimensionedScalar kAC
(
    parallelCompetitiveReaction.lookup("kAC")
);
