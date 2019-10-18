Info<< "Reading postProcessProperties\n" << endl;

IOdictionary postProcessProperties
(
    IOobject
    (
       "postProcessProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ,
       IOobject::NO_WRITE
    )
);

Info<< "Reading reactionProperties\n" << endl;

IOdictionary reactionType
(
    IOobject
    (
        "reactionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading reactionType\n" << endl;

const word reaction
(
    reactionType.lookup("reactionType")
);

vector center
(
    postProcessProperties.lookup("center")
);

vector centerAxis
(
    postProcessProperties.lookup("centerAxis")
);

dimensionedScalar bR
(
    postProcessProperties.lookup("bubbleRadius")
);

dimensionedScalar sectors
(
    postProcessProperties.lookup("sectors")
);

word interfacePatch
(
    postProcessProperties.lookup("interfacePatch")
);
