#include "initPhysisorption.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        #include "correctFluxesPhys.H"
        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phiA, A)
          - fvm::laplacian(A_Diff, A)
        );
        AEqn.solve();
    }
    runTime.write();
}
