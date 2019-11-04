// fill feature tensors and run forward pass
forAll(bouPatch, faceI)
{
    fTensorAB[faceI][0] = distF[faceI];
    fTensorAB[faceI][1] = A.internalField()[adjacentCellIDs[faceI]];
    fTensorAB[faceI][2] = B.internalField()[adjacentCellIDs[faceI]];
    fTensorAP[faceI][0] = distF[faceI];
    fTensorAP[faceI][1] = A.internalField()[adjacentCellIDs[faceI]];
    fTensorAP[faceI][2] = P.internalField()[adjacentCellIDs[faceI]];
    fTensorABP[faceI][0] = distF[faceI];
    fTensorABP[faceI][1] = A.internalField()[adjacentCellIDs[faceI]];
    fTensorABP[faceI][2] = B.internalField()[adjacentCellIDs[faceI]];
    fTensorABP[faceI][3] = P.internalField()[adjacentCellIDs[faceI]];
    fTensorABPS[faceI][0] = distF[faceI];
    fTensorABPS[faceI][1] = A.internalField()[adjacentCellIDs[faceI]];
    fTensorABPS[faceI][2] = B.internalField()[adjacentCellIDs[faceI]];
    fTensorABPS[faceI][3] = P.internalField()[adjacentCellIDs[faceI]];
    fTensorABPS[faceI][4] = S.internalField()[adjacentCellIDs[faceI]];
}

std::vector<torch::jit::IValue> IVal_fTensorABP{fTensorABP};
torch::Tensor lTensorA = A_model.forward(IVal_fTensorABP).toTensor();
auto lAccessorA = lTensorA.accessor<double,2>();
torch::Tensor lTensorB = B_model.forward(IVal_fTensorABP).toTensor();
auto lAccessorB = lTensorB.accessor<double,2>();
torch::Tensor lTensorP = P_model.forward(IVal_fTensorABP).toTensor();
auto lAccessorP = lTensorP.accessor<double,2>();
std::vector<torch::jit::IValue> IVal_fTensorABPS{fTensorABPS};
torch::Tensor lTensorS = S_model.forward(IVal_fTensorABPS).toTensor();
auto lAccessorS = lTensorS.accessor<double,2>();
std::vector<torch::jit::IValue> IVal_fTensorAB{fTensorAB};
torch::Tensor lTensorRAB = rAB_model.forward(IVal_fTensorAB).toTensor();
auto lAccessorRAB = lTensorRAB.accessor<double,1>();
std::vector<torch::jit::IValue> IVal_fTensorAP{fTensorAP};
torch::Tensor lTensorRAP = rAP_model.forward(IVal_fTensorAP).toTensor();
auto lAccessorRAP = lTensorRAP.accessor<double,1>();

// correct fluxes
const surfaceScalarField snGradA(fvc::snGrad(A));
const surfaceScalarField A_f(fvc::interpolate(A));
const surfaceScalarField snGradB(fvc::snGrad(B));
const surfaceScalarField B_f(fvc::interpolate(B));
const surfaceScalarField snGradP(fvc::snGrad(P));
const surfaceScalarField P_f(fvc::interpolate(P));
const surfaceScalarField snGradS(fvc::snGrad(S));
const surfaceScalarField S_f(fvc::interpolate(S));

scalarField sgsGradA_s(bouPatch.size(), 0.0);
scalarField linGradA_s(bouPatch.size(), 0.0);

forAll(bouPatch, faceI)
{
    // compute weights
    scalar w_num = Foam::pow(A.internalField()[adjacentCellIDs[faceI]], 2.0);
    scalar w_sgs = 1.0 - w_num;

    //
    // Correct species A
    //
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

    //
    // Correct species B
    //
    // correct gradient of B at first face opposite to boundary face
    scalar numGradB_f = snGradB[oppFaceIDs[faceI]];
    scalar ratioGradB_f = 1.0;
    if (mag(numGradB_f) > SMALL)
    {
        ratioGradB_f = (w_sgs * mag(lAccessorB[faceI][0]) + w_num * mag(numGradB_f)) / numGradB_f;
        //Info << "face I / corr: " << faceI << " / " << mag(ratioGradA_f) << "\n";
    }
    B_Diff.ref()[oppFaceIDs[faceI]] = DB.value() * mag(ratioGradB_f);

    // correct face value of B at first face opposite to boundary
    scalar numB_f = B_f[oppFaceIDs[faceI]];
    scalar sgsB_f = lAccessorB[faceI][1];
    sgsB_f = w_sgs * sgsB_f + w_num * numB_f;
    scalar B_f_min = max(min(B[adjacentCellIDs[faceI]], B[secAdjacentCellIDs[faceI]]), 0.0);
    scalar B_f_max = max(B[adjacentCellIDs[faceI]], B[secAdjacentCellIDs[faceI]]);
    bool correct_B_f =  (sgsB_f >= 0.0)
                     && (sgsB_f >= B_f_min)
                     && (sgsB_f <= B_f_max)
                     && (numB_f > SMALL);
    scalar ratioB_f = 1.0;
    if ( correct_B_f )
    {
        ratioB_f = sgsB_f / numB_f;
        //Info << "face I / corr: " << faceI << " / " << ratioA_f << "\n";
    }
    phiB.ref()[oppFaceIDs[faceI]] = phi[oppFaceIDs[faceI]] * ratioB_f;

    //
    // Correct species P
    //
    // correct gradient of P at first face opposite to boundary face
    scalar numGradP_f = snGradP[oppFaceIDs[faceI]];
    scalar ratioGradP_f = 1.0;
    if (mag(numGradP_f) > SMALL)
    {
        ratioGradP_f = (w_sgs * mag(lAccessorP[faceI][0]) + w_num * mag(numGradP_f)) / numGradP_f;
        //Info << "face I / corr: " << faceI << " / " << mag(ratioGradA_f) << "\n";
    }
    P_Diff.ref()[oppFaceIDs[faceI]] = DP.value() * mag(ratioGradP_f);

    // correct face value of P at first face opposite to boundary
    scalar numP_f = P_f[oppFaceIDs[faceI]];
    scalar sgsP_f = lAccessorP[faceI][1];
    sgsP_f = w_sgs * sgsP_f + w_num * numP_f;
    scalar P_f_min = max(min(P[adjacentCellIDs[faceI]], P[secAdjacentCellIDs[faceI]]), 0.0);
    scalar P_f_max = max(P[adjacentCellIDs[faceI]], P[secAdjacentCellIDs[faceI]]);
    bool correct_P_f =  (sgsP_f >= 0.0)
                     && (sgsP_f >= P_f_min)
                     && (sgsP_f <= P_f_max)
                     && (numP_f > SMALL);
    scalar ratioP_f = 1.0;
    if ( correct_P_f )
    {
        ratioP_f = sgsP_f / numP_f;
        //Info << "face I / corr: " << faceI << " / " << ratioA_f << "\n";
    }
    phiP.ref()[oppFaceIDs[faceI]] = phi[oppFaceIDs[faceI]] * ratioP_f;

    //
    // Correct species S
    //
    // correct gradient of S at first face opposite to boundary face
    scalar numGradS_f = snGradS[oppFaceIDs[faceI]];
    scalar ratioGradS_f = 1.0;
    if (mag(numGradS_f) > SMALL)
    {
        ratioGradS_f = (w_sgs * mag(lAccessorS[faceI][0]) + w_num * mag(numGradS_f)) / numGradS_f;
        //Info << "face I / corr: " << faceI << " / " << mag(ratioGradA_f) << "\n";
    }
    S_Diff.ref()[oppFaceIDs[faceI]] = DS.value() * mag(ratioGradS_f);

    // correct face value of P at first face opposite to boundary
    scalar numS_f = S_f[oppFaceIDs[faceI]];
    scalar sgsS_f = lAccessorS[faceI][1];
    sgsS_f = w_sgs * sgsS_f + w_num * numS_f;
    scalar S_f_min = max(min(S[adjacentCellIDs[faceI]], S[secAdjacentCellIDs[faceI]]), 0.0);
    scalar S_f_max = max(S[adjacentCellIDs[faceI]], S[secAdjacentCellIDs[faceI]]);
    bool correct_S_f =  (sgsS_f >= 0.0)
                     && (sgsS_f >= S_f_min)
                     && (sgsS_f <= S_f_max)
                     && (numS_f > SMALL);
    scalar ratioS_f = 1.0;
    if ( correct_S_f )
    {
        ratioS_f = sgsS_f / numS_f;
        //Info << "face I / corr: " << faceI << " / " << ratioA_f << "\n";
    }
    phiS.ref()[oppFaceIDs[faceI]] = phi[oppFaceIDs[faceI]] * ratioS_f;

    //
    // Correct reaction source term A+B->P
    //
    scalar AB_num = A[adjacentCellIDs[faceI]] * B[adjacentCellIDs[faceI]];
    scalar AB_sgs = lAccessorRAB[faceI];
    scalar ratioRAB = 1.0;
    if (AB_num > SMALL)
    {
        ratioRAB = (w_sgs * AB_sgs + w_num * AB_num) / AB_num;
    }
    rAB.ref()[adjacentCellIDs[faceI]] *= ratioRAB;

    //
    // Correct reaction source term A+P->S
    //
    scalar AP_num = A[adjacentCellIDs[faceI]] * P[adjacentCellIDs[faceI]];
    scalar AP_sgs = lAccessorRAP[faceI];
    scalar ratioRAP = 1.0;
    if (AP_num > SMALL)
    {
        ratioRAP = (w_sgs * AP_sgs + w_num * AP_num) / AP_num;
    }
    rAP.ref()[adjacentCellIDs[faceI]] *= ratioRAP;
}

Info << "Global Sh_A sgs/lin: " << (sum(sgsGradA_s) / sum(magSf))
     << " / " << (sum(linGradA_s) / sum(magSf)) << "\n";
