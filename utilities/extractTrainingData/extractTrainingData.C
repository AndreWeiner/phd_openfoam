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
    extractTrainingData

Description

\*---------------------------------------------------------------------------*/

#include <vector>
#include "fvCFD.H"

std::vector<scalar> compute_volume_average(
    std::vector<label> cells,
    const scalarField &cellVolumnes,
    const scalarField &field
){
    std::vector<scalar> average(cells.size()-1);
    scalar volSum = 0.0;
    scalar fieldVolSum = 0.0;
    for (size_t cellI=0; cellI < cells.size()-1; ++cellI)
    {
        volSum += cellVolumnes[cells[cellI]];
        fieldVolSum += cellVolumnes[cells[cellI]]
                     * field[cells[cellI]];
        average[cellI] = fieldVolSum / volSum;
    }
    return average;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "maxDist",
        "scalar|'maxDist'",
        "Consider only cell faces with a distance smaller than maxDist"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initReaction.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar maxDist = args.opt<scalar>("maxDist");
    Info << "Maximum distance set to " << maxDist << "\n";
    #include "findRelevantFaces.H"

    OFstream outputFile(runTime.path()/"trainingData.csv");
    outputFile.precision(15);
    std::string header = "# t, x, y, dist";
    for (int fI=0; fI < numOfFields; ++fI)
    {
        header += ", " + fieldNames[fI] + "_av"
               +  ", " + fieldNames[fI] + "_f"
               +  ", grad" + fieldNames[fI] + "_f";
        if (fI==0)
        {
            header += ", grad" + fieldNames[fI] + "_s";
        }

    }
    if (reaction == "singleReaction")
    {
        header += ", rAB_av";
    }
    else if(reaction == "parallelConsecutiveReaction")
    {
        header += ", rAB_av, rAP_av";
    }
    outputFile << header;

    instantList timeDirs = timeSelector::select0(runTime, args);

    for (int timei=1; timei<timeDirs.size(); ++timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        #include "createSpeciesFields.H"
/*
        volScalarField cellInd
        (
            IOobject
            (
              "cellInd",
              runTime.timeName(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("", dimensionSet(0,-3,0,1,0,0,0), 0.0)
        );

        forAll(mesh.C(), cellI)
        {
            cellInd.ref()[cellI] = 0.0;
        }
*/
        // loop over all fields, patch faces, and faces normal to patch face
        std::vector<std::vector<std::vector<scalar>>> C_av(
            numOfFields, std::vector<std::vector<scalar>>(faceIDs.size())
        );
        std::vector<std::vector<std::vector<scalar>>> gradC_f(
            numOfFields, std::vector<std::vector<scalar>>(faceIDs.size())
        );
        std::vector<std::vector<std::vector<scalar>>> C_f(
            numOfFields, std::vector<std::vector<scalar>>(faceIDs.size())
        );
        std::vector<std::vector<scalar>> gradC_s(
            numOfFields, std::vector<scalar>(faceIDs.size(), 0.0)
        );
        forAll (fields, fI)
        {
            const surfaceScalarField Fi(fvc::interpolate(fields[fI]));
            const surfaceVectorField gradFi(fvc::interpolate(fvc::grad(fields[fI])));
            for (size_t pFaceI=0; pFaceI<faceIDs.size(); ++pFaceI)
            {
                gradC_s[fI][pFaceI] = gradFi.boundaryField()[surfaceID][pFaceI]
                                    & surfaceNormal[pFaceI];
                std::vector<scalar> Ci_av = compute_volume_average(
                    cellIDs[pFaceI], mesh.V(), fields[fI].internalField()
                );

                std::vector<scalar> gradCi_f(faceIDs[pFaceI].size()-1, 0.0);
                std::vector<scalar> Ci_f(faceIDs[pFaceI].size()-1, 0.0);
                for (size_t bFaceI=0; bFaceI < faceIDs[pFaceI].size()-1; ++bFaceI)
                {
                    gradCi_f[bFaceI] = gradFi.internalField()[faceIDs[pFaceI][bFaceI+1]]
                                     & surfaceNormal[pFaceI];
                    Ci_f[bFaceI] = Fi.internalField()[faceIDs[pFaceI][bFaceI+1]];
                  //  cellInd.ref()[cellIDs[pFaceI][bFaceI]] = Ci_av[bFaceI];
                }
                C_av[fI][pFaceI] = Ci_av;
                gradC_f[fI][pFaceI] = gradCi_f;
                C_f[fI][pFaceI] = Ci_f;
            }
        }
        // write to file
        volScalarField AB("", fields[0]);
        volScalarField AP("", fields[0]);
        if (reaction == "singleReaction")
        {
            AB *= fields[1];
        }
        else if(reaction == "parallelConsecutiveReaction")
        {
            AB *= fields[1];
            AP *= fields[2];
        }

        for (size_t pFaceI=0; pFaceI < faceIDs.size(); ++pFaceI)
        {
            std::vector<label> bFaces = faceIDs[pFaceI];
            std::vector<scalar> AB_av;
            std::vector<scalar> AP_av;
            if (reaction == "singleReaction")
            {
                AB_av = compute_volume_average(
                    cellIDs[pFaceI], mesh.V(), AB
                );
            }
            else if(reaction == "parallelConsecutiveReaction")
            {
                AB_av = compute_volume_average(
                    cellIDs[pFaceI], mesh.V(), AB
                );
                AP_av = compute_volume_average(
                    cellIDs[pFaceI], mesh.V(), AP
                );
            }
            for (size_t bFaceI=0; bFaceI < bFaces.size()-1; ++bFaceI)
            {
                outputFile << "\n"
                    << runTime.timeName() << ", "
                    << surfaceCf[pFaceI].x() << ", "
                    << surfaceCf[pFaceI].y() << ", "
                    << distCf[pFaceI][bFaceI+1];
                forAll(fields, fI)
                {
                    outputFile << ", " << C_av[fI][pFaceI][bFaceI]
                               << ", " << C_f[fI][pFaceI][bFaceI]
                               << ", " << gradC_f[fI][pFaceI][bFaceI];
                               if(fI == 0)
                               {
                                  outputFile << ", " << gradC_s[fI][pFaceI];
                               }
                } // loop over all fields
                if (reaction == "singleReaction")
                {
                    outputFile << ", " << AB_av[bFaceI];
                }
                else if(reaction == "parallelConsecutiveReaction")
                {
                    outputFile << ", " << AB_av[bFaceI]
                               << ", " << AP_av[bFaceI];
                }
                // cellInd.ref()[cellIDs[pFaceI][bFaceI]] = AB_av[bFaceI];

            } // loop over all faces normal to boundary face
        } // loop over all boundary faces
       // cellInd.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
