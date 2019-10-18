/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    bubbleSurfaceFields

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    instantList timeDirs = timeSelector::select0(runTime, args);

    // determine patch ID
    label surfaceID(-1);
    forAll (mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].name() == "bubble")
        {
            surfaceID = patchI;
        }
    }

    const vectorField Cf(mesh.Cf().boundaryField()[surfaceID]);
    const vectorField Sf(mesh.Sf().boundaryField()[surfaceID]);
    const vectorField nf(-Sf / mag(Sf));

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        // read fields
        volVectorField U
        (
          IOobject
          (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
          ),
          mesh
        );
        #include "createPhi.H"
        U.correctBoundaryConditions();

        volScalarField p
        (
          IOobject
          (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
          ),
          mesh
        );
        p.correctBoundaryConditions();

        // values at the bubble surface
        const surfaceScalarField pf(fvc::interpolate(p));
        const surfaceTensorField gradUf = fvc::interpolate(fvc::grad(U));
        const vectorField Uf = U.boundaryField()[surfaceID];


        // write surface data to file
        OFstream outputFile(runTime.path()/runTime.timeName()/"surfaceFields.csv");
        outputFile.precision(12);
        outputFile << "# x, y, Af, U_tau, pf, d_tau_U_tau";

        forAll(nf, faceI)
        {
            outputFile << "\n";
            // tangential vector
            vector tau = vector(nf[faceI].y(), -nf[faceI].x(), nf[faceI].z());
            outputFile << Cf[faceI].x() << ", "        // 1
                       << Cf[faceI].y() << ", "        // 2
                       << mag(Sf[faceI]) << ", "       // 3
                       << (Uf[faceI] & tau) << ", "    // 4
                       << pf.boundaryField()[surfaceID][faceI] << ", " // 5
                       << ((gradUf.boundaryField()[surfaceID][faceI] & tau) & tau); // 6
        }

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
