#include "initPhysisorption.H"
while (simple.loop())
{
    Info<< "Time = " << runTime.timeName() << nl << endl;

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix AEqn
        (
            fvm::ddt(A)
          + fvm::div(phi, A)
          - fvm::laplacian(DA, A)
        );
        AEqn.solve();
    }

    runTime.write();
}
