#include "initDecayReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        volScalarField rA("rA", kA*A);

        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
         ==
            -rA
        );

        fvScalarMatrix PEqn
        (
            fvm::ddt(P)
          + fvm::div(phi, P)
          - fvm::laplacian(DP, P)
         ==
            rA
        );

        AEqn.solve();
        PEqn.solve();
    }

    runTime.write();
}
