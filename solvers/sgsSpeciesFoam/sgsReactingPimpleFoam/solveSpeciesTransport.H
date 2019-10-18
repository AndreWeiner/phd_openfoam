
    if (reaction == "noReaction")
    {
        volScalarField& A = fields[0];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A) + fvm::div(phi, A) - fvm::laplacian(diff[0], A)
            );
            AEqn.solve();
        }
    }
    else if(reaction == "decayReaction")
    {
        volScalarField rA("rA", kxx[0]*Foam::pow(fields[0],alpha[0].value()));

        volScalarField& A = fields[0];
        volScalarField& P = fields[1];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A) + fvm::div(phiSgs, A) - fvm::laplacian(DiffSgs, A) ==
                  - alpha[0]*rA
            );
            AEqn.solve();

            fvScalarMatrix PEqn
            (
                fvm::ddt(P) + fvm::div(phi, P) - fvm::laplacian(diff[1], P) ==
                  alpha[0]*rA
            );
            PEqn.solve();
        }
    }
    else if(reaction == "singleReaction")
    {
        volScalarField rAB("rAB", kxx[0]*
                           Foam::pow(fields[0], alpha[0].value())*
                           Foam::pow(fields[1], beta[0].value()));

        volScalarField& A = fields[0];
        volScalarField& B = fields[1];
        volScalarField& P = fields[2];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A) + fvm::div(phi, A) - fvm::laplacian(diff[0], A) ==
                  -alpha[0]*rAB
            );
            AEqn.solve();

            fvScalarMatrix BEqn
            (
                fvm::ddt(B) + fvm::div(phi, B) - fvm::laplacian(diff[1], B) ==
                  -beta[0]*rAB
            );
            BEqn.solve();

            fvScalarMatrix PEqn
            (
                fvm::ddt(P) + fvm::div(phi, P) - fvm::laplacian(diff[2], P) ==
                  rAB
            );
            PEqn.solve();
        }

    }
    else if(reaction == "consecutiveReaction")
    {
        volScalarField rAB("rAB", kxx[0]*
                           Foam::pow(fields[0], alpha[0].value())*
                           Foam::pow(fields[1], beta[0].value()));
        volScalarField rAP("rAP", kxx[1]*
                           Foam::pow(fields[0], alpha[1].value())*
                           Foam::pow(fields[2], beta[1].value()));

        volScalarField& A = fields[0];
        volScalarField& B = fields[1];
        volScalarField& P = fields[2];
        volScalarField& S = fields[3];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A) + fvm::div(phi, A) - fvm::laplacian(diff[0], A) ==
                  -alpha[0]*rAB-alpha[1]*rAP
            );
            AEqn.solve();

            fvScalarMatrix BEqn
            (
                fvm::ddt(B) + fvm::div(phi, B) - fvm::laplacian(diff[1], B) ==
                  -beta[0]*rAB
            );
            BEqn.solve();

            fvScalarMatrix PEqn
            (
                fvm::ddt(P) + fvm::div(phi, P) - fvm::laplacian(diff[2], P) ==
                  rAB-beta[1]*rAP
            );
            PEqn.solve();

            fvScalarMatrix SEqn
            (
                fvm::ddt(S) + fvm::div(phi, S) - fvm::laplacian(diff[3], S) ==
                  rAP
            );
            SEqn.solve();
        }

    }
    else if(reaction == "competitiveReaction")
    {
        volScalarField rAB("rAB", kxx[0]*
                           Foam::pow(fields[0], alpha[0].value())*
                           Foam::pow(fields[1], beta[0].value()));
        volScalarField rAC("rAC", kxx[1]*
                           Foam::pow(fields[0], alpha[1].value())*
                           Foam::pow(fields[2], beta[1].value()));

        volScalarField& A = fields[0];
        volScalarField& B = fields[1];
        volScalarField& C = fields[2];
        volScalarField& P = fields[3];
        volScalarField& S = fields[4];

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix AEqn
            (
                fvm::ddt(A) + fvm::div(phi, A) - fvm::laplacian(diff[0], A) ==
                  -alpha[0]*rAB-alpha[1]*rAC
            );
            AEqn.solve();

            fvScalarMatrix BEqn
            (
                fvm::ddt(B) + fvm::div(phi, B) - fvm::laplacian(diff[1], B) ==
                  -beta[0]*rAB
            );
            BEqn.solve();

            fvScalarMatrix CEqn
            (
                fvm::ddt(C) + fvm::div(phi, C) - fvm::laplacian(diff[2], C) ==
                  -beta[1]*rAC
            );
            CEqn.solve();

            fvScalarMatrix PEqn
            (
                fvm::ddt(P) + fvm::div(phi, P) - fvm::laplacian(diff[3], P) == rAB
            );
            PEqn.solve();

            fvScalarMatrix SEqn
            (
                fvm::ddt(S) + fvm::div(phi, S) - fvm::laplacian(diff[4], S) == rAC
            );
            SEqn.solve();
        }

    }

    forAll(fields, i)
    {
      Foam::bound(fields[i], minC);
    }
