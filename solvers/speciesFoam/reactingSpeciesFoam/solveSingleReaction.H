#include "initSingleReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {

        volScalarField rAB("rAB",kAB*A*B);

        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
         ==
            -rAB
        );

        fvScalarMatrix AphyEqn
        (
            fvm::ddt(Aphy)
          + fvm::div(phi, Aphy)
          - fvm::laplacian(DA, Aphy)
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
            rAB
        );

        AEqn.solve();
        AphyEqn.solve();
        BEqn.solve();
        PEqn.solve();
    }

    runTime.write();
}
