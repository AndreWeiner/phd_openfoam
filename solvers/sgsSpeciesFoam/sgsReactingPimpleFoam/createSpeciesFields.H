    IOdictionary reactionProperties
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
        reactionProperties.lookup("reactionType")
    );

    int numOfFields = 0;
    List<word> fieldNames(5);
    std::vector<dimensionedScalar> diff{};
    std::vector<dimensionedScalar> kxx{};
    std::vector<dimensionedScalar> alpha{};
    std::vector<dimensionedScalar> beta{};

    if (reaction == "noReaction")
    {
        numOfFields = 1;
        fieldNames[0] = "A";
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DA")));
    }
    else if(reaction == "decayReaction")
    {
        numOfFields = 2;
        fieldNames[0] = "A";
        fieldNames[1] = "P";
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DA")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DP")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kA")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_1")));
    }
    else if(reaction == "singleReaction")
    {
        numOfFields = 3;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "P";
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DA")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DB")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DP")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kAB")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_1")));
        beta.push_back(dimensionedScalar(reactionProperties.lookup("beta_1")));
    }
    else if(reaction == "consecutiveReaction")
    {
        numOfFields = 4;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "P";
        fieldNames[3] = "S";
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DA")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DB")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DP")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DS")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kAB")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kAP")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_1")));
        beta.push_back(dimensionedScalar(reactionProperties.lookup("beta_1")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_2")));
        beta.push_back(dimensionedScalar(reactionProperties.lookup("beta_2")));
    }
    else if(reaction == "competitiveReaction")
    {
        numOfFields = 5;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "C";
        fieldNames[3] = "P";
        fieldNames[4] = "S";
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DA")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DB")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DC")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DP")));
        diff.push_back(dimensionedScalar(reactionProperties.lookup("DS")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kAB")));
        kxx.push_back(dimensionedScalar(reactionProperties.lookup("kAC")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_1")));
        beta.push_back(dimensionedScalar(reactionProperties.lookup("beta_1")));
        alpha.push_back(dimensionedScalar(reactionProperties.lookup("alpha_2")));
        beta.push_back(dimensionedScalar(reactionProperties.lookup("beta_2")));
    }
    else
    {
      Info << "Unknown reaction type " << reaction << endl;
      return 0;
    }

    PtrList<volScalarField> fields(numOfFields);

    forAll(fields, i)
    {
       fields.set
       (
           i,
           new volScalarField
           (
               IOobject
               (
                  fieldNames[i],
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE
               ),
           mesh
           )
       );
    }

    dimensionedScalar minC ("minC", fields[0].dimensions(), 0.0);

    Info << "Finished reading reaction properties for " << reaction << endl;
