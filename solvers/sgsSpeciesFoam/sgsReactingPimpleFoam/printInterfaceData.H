
vectorField cellC = mesh.C().boundaryField()[surfaceID];

scalarField CiBoundary = fields[0].boundaryField()[surfaceID];
scalarField localSh(CiBoundary.size(), 0.0);
scalar globalSh = 0.0;

// Linear gradient in case the SGS model is not applied or it does not converge
//scalarField nGradCi(fields[0].boundaryField()[surfaceID].patch().size(), 0.0);
//nGradCi = fields[0].boundaryField()[surfaceID].snGrad();
localSh = dcSigFinal / CiBoundary*refLength.value();

globalSh = sum(localSh*mag(mesh.Sf().boundaryField()[surfaceID]))
          /sum(mag(mesh.Sf().boundaryField()[surfaceID]));
Info << "Global Sherwood number: " << mag(globalSh) << endl;

// Get the polar angles
scalarField radii = mag(cellC);
scalarField polarAngles(cellC.size());
// Rising in Y-direction:
forAll (radii, i) {
  polarAngles[i] = Foam::acos(cellC[i].y()/radii[i])
                 * 180.0/pi
                 * (cellC[i].x()/mag(cellC[i].x()));
}

//extension for reaction:
//print interface concentration for each species
//and local Sherwood number for transfer species
OFstream of(runTime.timeName()/"interfaceData.dat");
Info << of << endl;
of.precision(12);

of  << "Time: " << runTime.timeName() << tab << "Global Sh: " << globalSh << endl
    << "# polarAngle" << tab << "dcSig" << tab << "dcBulk" << tab << "cBulk" << tab
    << "corrAdv" << tab << "corrFace" << tab << "corrSig" << endl;
/*
if (reaction == "noReaction")
    {
      of << "# y / Sh_A / c_A" << endl;
    }
    else if(reaction == "decayReaction")
    {
      of << "# y / Sh_A / c_A / c_P" << endl;
    }
    else if(reaction == "singleReaction")
    {
      of << "# y / Sh_A / c_A / c_B / c_P" << endl;
    }
    else if(reaction == "consecutiveReaction")
    {
      of << "# y / Sh_A / c_A / c_B / c_P / c_S" << endl;
    }
    else if(reaction == "competitiveReaction")
    {
      of << "# y / Sh_A / c_A / c_B / c_C / c_P / c_S" << endl;
    }

//save interface concentration for all species
List<scalarField> interfaceConc(numOfFields);

forAll(fields, i)
{
  interfaceConc[i] = fields[i].boundaryField()[surfaceID];
}
*/
//print
forAll(cellC, i)
{
    scalar angle = polarAngles[i];
    scalar locDcSig = localSh[i];
    scalar locDcBulk = dcBulkFinal[i];
    scalar locCBulk = cBulkFinal[i];
    scalar locCorrAdv = corrAdv[i];
    scalar locCorrFace = corrFace[i];
    scalar locCorrSig = corrSig[i];

    of  << angle << tab << locDcSig << tab << locDcBulk << tab << locCBulk << tab
    << locCorrAdv << tab << locCorrFace << tab << locCorrSig << endl;

/*    forAll(fields, k)
    {
	of << interfaceConc[k][i];

    }

    of << endl;*/
}
