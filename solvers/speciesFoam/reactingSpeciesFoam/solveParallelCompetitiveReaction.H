#include "initParallelCompetitiveReaction.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {

        volScalarField rAB("rAB",kAB*A*B);
        volScalarField rAC("rAC",kAC*A*C);

        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
         ==
            -(rAB+rAC)
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

        fvScalarMatrix CEqn
        (
            fvm::ddt(C)
          + fvm::div(phi, C)
          - fvm::laplacian(DC, C)
         ==
            -rAC
        );

        fvScalarMatrix PEqn
        (
            fvm::ddt(P)
          + fvm::div(phi, P)
          - fvm::laplacian(DP, P)
         ==
            rAB
        );

        fvScalarMatrix SEqn
        (
            fvm::ddt(S)
          + fvm::div(phi, S)
          - fvm::laplacian(DS, S)
         ==
            rAC
        );

        AEqn.solve();
        AphyEqn.solve();
        BEqn.solve();
        CEqn.solve();
        PEqn.solve();
        SEqn.solve();
    }

    runTime.write();
}
