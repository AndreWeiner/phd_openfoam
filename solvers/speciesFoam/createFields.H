Info<< "Reading surface flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading reactionType\n" << endl;

IOdictionary reactionType
(
    IOobject
    (
        "reactionType",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word reaction
(
    reactionType.lookup("reactionType")
);
