Info << "Parallel consecutive reaction:" << endl;
Info << "A + B -(kAB)> P" << endl;
Info << "P + B -(kPB)> S" << endl;
Info << "A - Transfer species" << endl;
Info << "B - Bulk species" << endl;
Info << "Reading transport and reaction properties." << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
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
