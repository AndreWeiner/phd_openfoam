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

int numOfFields = 0;
List<word> fieldNames(10);

if (reaction == "noReaction")
{
    numOfFields = 1;
    fieldNames[0] = "A";
}
else if(reaction == "decayReaction")
{
    numOfFields = 2;
    fieldNames[0] = "A";
    fieldNames[1] = "P";
}
else if(reaction == "singleReaction")
{
    numOfFields = 3;
    fieldNames[0] = "A";
    fieldNames[1] = "B";
    fieldNames[2] = "P";
}
else if(reaction == "consecutiveReaction")
{
    numOfFields = 4;
    fieldNames[0] = "A";
    fieldNames[1] = "B";
    fieldNames[2] = "P";
    fieldNames[3] = "S";
}
else
{
  Info << "Unknown reaction type " << reaction << endl;
  return 0;
}
