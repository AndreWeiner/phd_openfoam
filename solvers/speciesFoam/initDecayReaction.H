Info << "Decay reaction:" << endl;
Info << "A -(kA)> P" << endl;
Info << "A - Transfer species" << endl;
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
