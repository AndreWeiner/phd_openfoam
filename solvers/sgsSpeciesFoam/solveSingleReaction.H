#include "initSingleReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        volScalarField rAB("rAB", kAB*A*B);

        #include "correctFluxesSingle.H"
        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phiA, A)
          - fvm::laplacian(A_Diff, A)
         ==
            -rAB
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
            rAB
        );

        AEqn.solve();
        BEqn.solve();
        PEqn.solve();
    }

    runTime.write();
}
