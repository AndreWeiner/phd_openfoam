#include "initParallelConsecutiveReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        volScalarField rAB("rAB",kAB*A*B);
        volScalarField rAP("rAP",kAP*A*P);

        #include "correctFluxesConsecutive.H"
        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phiA, A)
          - fvm::laplacian(A_Diff, A)
         ==
            -(rAB+rAP)
        );

        fvScalarMatrix BEqn
        (
            fvm::ddt(B)
          + fvm::div(phiB, B)
          - fvm::laplacian(B_Diff, B)
         ==
            -rAB
        );

        fvScalarMatrix PEqn
        (
            fvm::ddt(P)
          + fvm::div(phiP, P)
          - fvm::laplacian(P_Diff, P)
         ==
            rAB-rAP
        );

        fvScalarMatrix SEqn
        (
            fvm::ddt(S)
          + fvm::div(phiS, S)
          - fvm::laplacian(S_Diff, S)
         ==
            rAP
        );

        AEqn.solve();
        BEqn.solve();
        PEqn.solve();
        SEqn.solve();
    }

    runTime.write();
}
