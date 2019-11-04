#include "initDecayReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {

        volScalarField rA("rA",kA*A);

        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
         ==
            -rA
        );

        fvScalarMatrix AphyEqn
        (
            fvm::ddt(Aphy)
          + fvm::div(phi, Aphy)
          - fvm::laplacian(DA, Aphy)
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
        AphyEqn.solve();
        PEqn.solve();
    }

    runTime.write();
}
