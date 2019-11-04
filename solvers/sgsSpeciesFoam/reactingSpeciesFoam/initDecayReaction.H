Info << "Decay reaction:" << endl;
Info << "A -(kA)> P" << endl;
Info << "A - Transfer species" << endl;
Info << "Reading decayReaction properties." << endl;

IOdictionary decayReaction
(
    IOobject
    (
        "decayReaction",
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
    decayReaction.lookup("DA")
);

Info<< "Reading diffusivity DP\n" << endl;

dimensionedScalar DP
(
    decayReaction.lookup("DP")
);

Info<< "Reading reaction rate kA\n" << endl;

dimensionedScalar kA
(
    decayReaction.lookup("kA")
);
