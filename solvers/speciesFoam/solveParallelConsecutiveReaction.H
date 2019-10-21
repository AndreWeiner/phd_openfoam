#include "initParallelConsecutiveReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        volScalarField rAB("rAB",kAB*A*B);
        volScalarField rAP("rAP",kAP*A*P);

        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
         ==
            -(rAB+rAP)
        );

        fvScalarMatrix BEqn
        (
            fvm::ddt(B)
          + fvm::div(phi, B)
          - fvm::laplacian(DB, B)
         ==
            -rAB
        );

        fvScalarMatrix PEqn
        (
            fvm::ddt(P)
          + fvm::div(phi, P)
          - fvm::laplacian(DP, P)
         ==
            rAB-rAP
        );

        fvScalarMatrix SEqn
        (
            fvm::ddt(S)
          + fvm::div(phi, S)
          - fvm::laplacian(DS, S)
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
