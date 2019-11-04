#include "initDecayReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        volScalarField rA("rA", kA*A);

        #include "correctFluxesDecay.H"
        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phiA, A)
          - fvm::laplacian(A_Diff, A)
         ==
            -rA
        );

        fvScalarMatrix PEqn
        (
            fvm::ddt(P)
          + fvm::div(phiP, P)
          - fvm::laplacian(P_Diff, P)
         ==
            rA
        );

        AEqn.solve();
        PEqn.solve();
    }

    runTime.write();
}
