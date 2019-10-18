// Vol field: bulk diffusivity, to be interpolated to the faces
volScalarField DiffSgsVF
(
    IOobject
    (
        "DiffSgsVF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    diff[0]
);

// Surface field: interpolated diffusivities from cell centre to face centre
// To be corrected from the SGS model
surfaceScalarField DiffSgs
(
    IOobject
    (
        "DiffSgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(DiffSgsVF),
    fixedValueFvPatchScalarField::typeName
);

// Flux to be corrected from the SGS model
surfaceScalarField phiSgs
(
    IOobject
    (
        "phiSgs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phi
);

// gradient comming from standard interpolation
surfaceScalarField snGradNum
(
    IOobject
    (
        "snGradNum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::snGrad(fields[0])
);

// face value comming from standard interpolation
surfaceScalarField cfNum
(
    IOobject
    (
        "cfNum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(fields[0])
);

// save face distance for later use
scalarField rf(bouPatch.size(),0.0);

if (applySGS)
{
  // fill feature tensor
  int nCellsExtrapolation = 0;
  //bool extrapolate = false;
  scalarField extrapolate(bouPatch.size(),0.0);
  forAll (bouPatch, faceI)
  {
      rf[faceI] = mag(mesh.faceCentres()[patchFaceIDs[faceI]] - mesh.faceCentres()[oppFaceIDs[faceI]]);
      scalar cNorm = (fields[0][adjacentCellIDs[faceI]] - featureMinValues[0])
                       / (featureMaxValues[0] - featureMinValues[0]);
      featureTensor[faceI][1] = (fields[0][adjacentCellIDs[faceI]] - featureMinValues[0])
                              / (featureMaxValues[0] - featureMinValues[0]);
      featureTensor[faceI][0] = (rf[faceI] - featureMinValues[1])
                              / (featureMaxValues[1] - featureMinValues[1]);
      featureTensor[faceI][2] = (kxx[0].value() * refLength.value() / riseVelocity.value() - featureMinValues[2])
                              / (featureMaxValues[2] - featureMinValues[2]);
      if (cNorm < SMALL)
      {
        extrapolate[faceI] = 1.0;
        nCellsExtrapolation++;
      }
  }

  if (nCellsExtrapolation > 0)
  {
    std::cout << "Warning: extrapolating feature space for " << nCellsExtrapolation
    << " out of " << nInterfaceCells << " cells." << std::endl;
  }

    // convert tensor to device-specific input value
    std::vector<torch::jit::IValue> modelFeatures{featureTensor};

    // compute and correct bulk face gradient
    if (diffFacesCorr)
    {
        torch::Tensor labelTwoTensor = dc_bulk_module->forward(modelFeatures).toTensor();
        auto labelTwoAccessor = labelTwoTensor.accessor<float,2>();

        forAll(oppFaceIDs, faceI)
        {
          // convert label to actual concentration gradient
          double faceGradient = fields[0][adjacentCellIDs[faceI]] / rf[faceI] * (
                  labelMinValues[1] + labelTwoAccessor[faceI][0] * (labelMaxValues[1] - labelMinValues[1])
          );

          //exception handling if gradient is too small
          if ( (mag(snGradNum[oppFaceIDs[faceI]]) > 1.0) && (mag(faceGradient) > SMALL) )
          {
            if (blendDiffBulk) //blending
            {
              faceGradient = (1.0 - Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaDiff.value())) * mag(faceGradient)
                           + Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaDiff.value()) * mag(snGradNum[oppFaceIDs[faceI]]);
            }

            dcBulkFinal[faceI] = faceGradient;
            scalar corrFact = mag(faceGradient
                     / snGradNum[oppFaceIDs[faceI]]);
            corrFace[faceI]= corrFact;
            DiffSgs[oppFaceIDs[faceI]] = diff[0].value()*corrFact;
          }
          else
            {corrFace[faceI] =-1.0;}
        }
    }

    // compute and correct bulk face values
    if (advCorr)
    {
        torch::Tensor labelOneTensor = c_bulk_module->forward(modelFeatures).toTensor();
        auto labelOneAccessor = labelOneTensor.accessor<float,2>();

        forAll(oppFaceIDs, faceI)
        {
          // convert label to actual concentration value
          double faceValue = fields[0][adjacentCellIDs[faceI]] - (
                    labelMinValues[0] + labelOneAccessor[faceI][0] * (labelMaxValues[0] - labelMinValues[0])
          );

          scalar cfMin = fields[0][secAdjacentCellIDs[faceI]];
          scalar cfMax = fields[0][adjacentCellIDs[faceI]];
          bool inbetween = (faceValue >= cfMin) && (faceValue <= cfMax);
          if (!inbetween) //set concentration value to minimum
            {faceValue = cfMin;}
          //exception handling if concentration value is too small
          if ( (faceValue > SMALL) && (cfNum[oppFaceIDs[faceI]] > SMALL) )
          {
            if (blendAdvection) //blending
            {
              faceValue = (1.0 - Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaAdv.value())) * faceValue
                        + Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaAdv.value()) * cfNum[oppFaceIDs[faceI]];
            }

            cBulkFinal[faceI] = faceValue;
            scalar corrFact= faceValue / cfNum[oppFaceIDs[faceI]];
            corrAdv[faceI] = corrFact;
            phiSgs[oppFaceIDs[faceI]] = phi[oppFaceIDs[faceI]]
                                      * faceValue / cfNum[oppFaceIDs[faceI]];
          }
          else
            {corrAdv[faceI] =-1.0;}
        }
    }

    // compute and correct interface gradient
    if (diffSigmaCorr)
    {
        torch::Tensor labelThreeTensor = dc_sig_module->forward(modelFeatures).toTensor();
        auto labelThreeAccessor = labelThreeTensor.accessor<float,2>();

        forAll(bouPatch, faceI)
        {
            // convert label to actual concentration gradient
            double faceGradient = 1.0 / (fields[0][adjacentCellIDs[faceI]] * rf[faceI])
                  * (labelMinValues[2] + labelThreeAccessor[faceI][0] * (labelMaxValues[2] - labelMinValues[2]) );
            //exception handling if gradient is too small
            if ( (mag(snGradNum.boundaryField()[surfaceID][faceI]) > SMALL) && (mag(faceGradient) > SMALL) )
            {
              if (blendDiffSig) //blending
              {
                faceGradient = (1.0 - Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaDiffSigma.value())) * mag(faceGradient)
                             + Foam::pow(fields[0][adjacentCellIDs[faceI]],alphaDiffSigma.value())
                                      * mag(snGradNum.boundaryField()[surfaceID][faceI]);
              }
              dcSigFinal[faceI] = faceGradient;  //write final faceGradient to plot local Sh
              scalar corrFact = mag(faceGradient / snGradNum.boundaryField()[surfaceID][faceI]);
              corrSig[faceI] = corrFact;
              DiffSgs.boundaryFieldRef()[surfaceID][faceI] = diff[0].value()*corrFact;
            }
            else
            {
              corrSig[faceI] = -1.0;
              dcSigFinal[faceI] = -1.0;
            }
        }
    }

}
