// fill feature tensors and run forward pass
forAll(bouPatch, faceI)
{
    fTensorA[faceI][0] = distF[faceI];
    fTensorA[faceI][1] = A.internalField()[adjacentCellIDs[faceI]];
}

std::vector<torch::jit::IValue> IVal_fTensorA{fTensorA};
torch::Tensor lTensorA = A_model.forward(IVal_fTensorA).toTensor();
auto lAccessorA = lTensorA.accessor<double,2>();

// correct fluxes
const surfaceScalarField snGradA(fvc::snGrad(A));
const surfaceScalarField A_f(fvc::interpolate(A));

scalarField sgsGradA_s(bouPatch.size(), 0.0);
scalarField linGradA_s(bouPatch.size(), 0.0);

forAll(bouPatch, faceI)
{
    // compute weights
    scalar w_num = Foam::pow(A.internalField()[adjacentCellIDs[faceI]], 2.0);
    scalar w_sgs = 1.0 - w_num;
    // correct interface gradient of A
    scalar numGradA_s = snGradA.boundaryField()[surfaceID][faceI];
    scalar ratioGradA_s = 1.0;
    if (mag(numGradA_s) > SMALL)
    {
        ratioGradA_s = (w_sgs * mag(lAccessorA[faceI][0]) + w_num * mag(numGradA_s)) / numGradA_s;
        //ratioGradA_s = lAccessorA[faceI][0] / numGradA_s;
    }
    A_Diff.boundaryFieldRef()[surfaceID][faceI] = DA.value() * mag(ratioGradA_s);
    sgsGradA_s[faceI] = numGradA_s * ratioGradA_s * magSf[faceI];
    linGradA_s[faceI] = numGradA_s * magSf[faceI];

    // correct gradient of A at first face opposite to boundary face
    scalar numGradA_f = snGradA[oppFaceIDs[faceI]];
    scalar ratioGradA_f = 1.0;
    if (mag(numGradA_f) > SMALL)
    {
        ratioGradA_f = (w_sgs * mag(lAccessorA[faceI][1]) + w_num * mag(numGradA_f)) / numGradA_f;
        //ratioGradA_f = lAccessorA[faceI][1] / numGradA_f;
        //Info << "face I / corr: " << faceI << " / " << mag(ratioGradA_f) << "\n";
    }
    A_Diff.ref()[oppFaceIDs[faceI]] = DA.value() * mag(ratioGradA_f);

    // correct face value of A at first face opposite to boundary
    scalar numA_f = A_f[oppFaceIDs[faceI]];
    scalar sgsA_f = lAccessorA[faceI][2];
    sgsA_f = w_sgs * sgsA_f + w_num * numA_f;
    scalar A_f_min = max(min(A[adjacentCellIDs[faceI]], A[secAdjacentCellIDs[faceI]]), 0.0);
    scalar A_f_max = max(A[adjacentCellIDs[faceI]], A[secAdjacentCellIDs[faceI]]);
    bool correct_A_f =  (sgsA_f >= 0.0)
                     && (sgsA_f >= A_f_min)
                     && (sgsA_f <= A_f_max)
                     && (numA_f > SMALL);
    scalar ratioA_f = 1.0;
    if ( correct_A_f )
    {
        ratioA_f = sgsA_f / numA_f;
        //Info << "face I / corr: " << faceI << " / " << ratioA_f << "\n";
    }
    phiA.ref()[oppFaceIDs[faceI]] = phi[oppFaceIDs[faceI]] * ratioA_f;
}

Info << "Global Sh_A sgs/lin: " << (sum(sgsGradA_s) / sum(magSf))
     << " / " << (sum(linGradA_s) / sum(magSf)) << "\n";
