Info << "Physisorption:" << endl;
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

Info<< "Reading diffusivity DA\n" << endl;

dimensionedScalar DA
(
    transportProperties.lookup("DA")
);
