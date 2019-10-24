/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    localeSh

Description
    Compute species gradient and value at the interface.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initReaction.H"


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
    scalarField globalGrad(timeDirs.size(), 0.0);

    forAll (timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        List<scalarField> interfaceConc(numOfFields);

        forAll (fields, i)
        {
           interfaceConc[i] = fields[i].boundaryField()[surfaceID];
        }

        scalarField gradA(fields[0].boundaryField()[surfaceID].snGrad());
        globalGrad[timei] = sum(gradA * mag(Sf)) / sum(mag(Sf));

        OFstream outputFile(runTime.path()/runTime.timeName()/"interfaceData.csv");
        outputFile.precision(15);

        if (reaction == "noReaction")
        {
            outputFile  << "# x, y, A, snGradA";
        }
        else if(reaction == "decayReaction")
        {
            outputFile  << "# x, y, A, snGradA, c_P";
        }
        else if(reaction == "singleReaction")
        {
            outputFile  << "# x, y, A, snGradA, c_P, c_B";
        }
        else if(reaction == "parallelConsecutiveReaction")
        {
            outputFile  << "# x, y, A, snGradA, c_P, c_B, c_S";
        }

        forAll(Cf, faceI)
        {
            outputFile << "\n"
                      << Cf[faceI].x() << ", " << Cf[faceI].y() << ", "
                      << mag(Sf[faceI]) << ", " << gradA[faceI];
            forAll (fields, i)
            {
                if (i > 0)
                {
                    outputFile << ", " << fields[i][faceI];
                }
            }
        }
    } // end of time loop

    // write global gradient
    OFstream gradFile(runTime.path()/"gradGlobal.csv");
    gradFile.precision(15);

    gradFile << "# time, gradA_gl, A_gl";

    forAll (timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        gradFile << "\n"
                 << runTime.timeName() << ", "
                 << globalGrad[timei] << ", "
                 << sum(mag(Sf));
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
