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
    Compute the locale sherwood number

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    timeSelector::addOptions(true, false);
    argList args(argc, argv);
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"
    #include "createFields.H"

    label patchID = mesh.boundaryMesh().findPatchID(interfacePatch);

    vectorField faceCenters = mesh.Cf().boundaryField()[patchID];
    scalarField polarAngles(faceCenters.size());
    scalarField polarRadius(faceCenters.size());
    scalarField radii(faceCenters.size());
    scalarField cellFaces = mesh.magSf().boundaryField()[patchID];
    vector v1 = vector(0, 0, 0);
    vector v2 = vector(0, 0, 0);
    vector v3 = vector(0, 0, 0);
    if(centerAxis.x()>0)
    {
        v1 = vector(1, 0, 0);
        v2 = vector(0, 1, 0);
        v3 = vector(0, 0, 1);
    }
    else if(centerAxis.y()>0)
    {
        v1 = vector(0, 1, 0);
        v2 = vector(0, 0, 1);
        v3 = vector(1, 0, 0);
    }
    else if(centerAxis.z()>0)
    {
        v1 = vector(0, 0, 1);
        v2 = vector(1, 0, 0);
        v3 = vector(0, 1, 0);
    }


    forAll(faceCenters, faceI)
    {
        vector c = faceCenters[faceI] - center;
        scalar x = c & centerAxis;
        scalar y = c & v2;
        scalar z = c & v3;
        scalar yz = ::sqrt(::pow(y, 2.0) + ::pow(z, 2.0));
        scalar r = ::mag(c);
        radii[faceI] = yz;
        scalar t = ::asin(yz/r);
        if (x < 0.0)
        {
            t=constant::mathematical::pi-t;
        }
        polarAngles[faceI] = t;
        polarRadius[faceI] = yz*constant::mathematical::pi*2.0;
    }

    int numOfFields = 0;
    List<word> fieldNames(10);

    if (reaction == "noReaction")
    {
        numOfFields = 1;
        fieldNames[0] = "A";
    }
    else if(reaction == "decayReaction")
    {
        numOfFields = 2;
        fieldNames[0] = "A";
        fieldNames[1] = "P";
    }
    else if(reaction == "singleReaction")
    {
        numOfFields = 3;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "P";
    }
    else if(reaction == "consecutiveReaction")
    {
        numOfFields = 4;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "P";
        fieldNames[3] = "S";
    }
    else if(reaction == "competitiveReaction")
    {
        numOfFields = 5;
        fieldNames[0] = "A";
        fieldNames[1] = "B";
        fieldNames[2] = "C";
        fieldNames[3] = "P";
        fieldNames[4] = "S";
    }
    else
    {
      Info << "Unknown reaction type " << reaction << endl;
      return 0;
    }

    PtrList<volScalarField> fields(numOfFields);

    forAll(fields, i)
    {
       fields.set
       (
           i,
           new volScalarField
           (
               IOobject
               (
                  fieldNames[i],
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
               ),
           mesh
           )
       );
    }

    List<scalarField> interfaceConc(numOfFields);
    List<scalar> averageConc(numOfFields);

    forAll(fields, i)
    {
       interfaceConc[i] = fields[i].boundaryField()[patchID];
       averageConc[i] = sum(interfaceConc[i]*cellFaces)/sum(cellFaces);
    }

    scalarField gradA(fields[0].boundaryField()[patchID].snGrad());
    scalarField localShA(gradA/interfaceConc[0]*2.0*bR.value());
    scalar globalShA(sum(cellFaces*localShA)/sum(cellFaces));

    List<scalarField> sectorValues(numOfFields);
    scalarField sectorAngle(int(sectors.value()), 0.0);
    scalarField sectorRadius(int(sectors.value()), 0.0);
    scalarField sectorArea(int(sectors.value()), 0.0);
    scalarField sectorShA(int(sectors.value()), 0.0);
    scalarField secStart(int(sectors.value()), 0.0);
    scalarField secEnd(int(sectors.value()), 0.0);

    forAll(sectorValues, fieldI)
    {
      sectorValues[fieldI] = secStart*0.0;
    }

    forAll(secEnd, secI)
    {
      secStart[secI] = constant::mathematical::pi/sectors.value()*(secI);
      secEnd[secI] = constant::mathematical::pi/sectors.value()*(secI+1);
    }

    for (int i=0; i < sectors.value(); i++)
    {
      forAll(polarAngles, faceI)
      {
        if (polarAngles[faceI] >= secStart[i] && polarAngles[faceI] < secEnd[i])
        {
          forAll(interfaceConc, fieldI)
          {
            sectorValues[fieldI][i] += interfaceConc[fieldI][faceI]*cellFaces[faceI];
          } // loop over all concentrations
          sectorArea[i] += cellFaces[faceI];
          sectorAngle[i] += polarAngles[faceI]*cellFaces[faceI];
          sectorRadius[i] += polarRadius[faceI]*cellFaces[faceI];
          sectorShA[i] += localShA[faceI]*cellFaces[faceI];
        } // in sector?
      } // loop over faces
      if (sectorArea[i] > 0)
      {
        sectorAngle[i] /= sectorArea[i];
        sectorRadius[i] /= sectorArea[i];
        sectorShA[i] /= sectorArea[i];
        forAll(sectorValues, fieldI)
        {
          sectorValues[fieldI][i] /= sectorArea[i];
        }
      } // normalize with area if at least one cell is in the sector
    } // loop over sectors

    OFstream outputFile(runTime.path()/runTime.timeName()/"interfaceData.dat");
    outputFile.precision(12);

    outputFile  << "# time: " << runTime.timeName() << endl
                << "# effective area: " << sum(cellFaces) << endl
                << "# global Sh: " << tab << globalShA << endl;

    if (reaction == "noReaction")
    {
        outputFile << "# average c_A: " << averageConc[0] << endl
                   << "# polar Angle / axis Radius / Sh_A / c_A" << endl;
    }
    else if(reaction == "decayReaction")
    {
      outputFile << "# average c_A / c_P: "
                 << averageConc[0] << tab << averageConc[1] << endl
                 << "# polar Angle / axis Radius / Sh_A / c_A / c_P" << endl;
    }
    else if(reaction == "singleReaction")
    {
      outputFile << "# average c_A / c_B / c_P: "
                 << averageConc[0] << tab << averageConc[1] << tab
                 << averageConc[2] << endl
                 << "# polar Angle / axis Radius / Sh_A / c_A / c_B / c_P" << endl;
    }
    else if(reaction == "consecutiveReaction")
    {
      outputFile << "# average c_A / c_B / c_P / c_S: "
                 << averageConc[0] << tab << averageConc[1] << tab
                 << averageConc[2] << tab << averageConc[3] << endl
                 << "# polar Angle / axis Radius / Sh_A / c_A / c_B / c_P / c_S" << endl;
    }
    else if(reaction == "competitiveReaction")
    {
      outputFile << "# average c_A / c_B / c_C / c_P / c_S: "
                 << averageConc[0] << tab << averageConc[1] << tab
                 << averageConc[2] << tab << averageConc[3] << tab
                 << averageConc[4] << endl
                 << "# polar Angle / axis Radius / Sh_A / c_A / c_B / c_C / c_P / c_S" << endl;
    }

    forAll(sectorAngle, angleI)
    {
      outputFile  << sectorAngle[angleI] << tab << sectorRadius[angleI] << tab
                  << sectorShA[angleI] << tab;
      forAll(sectorValues, fieldI)
      {
        outputFile << sectorValues[fieldI][angleI] << tab;
      }
      outputFile << endl;
    }
/*    forAll(faceCenters,faceI)
    {
      outputFile  << polarAngles[faceI] << tab << polarRadius[faceI] << tab
                  << localShA[faceI] << tab;
      forAll(interfaceConc, i)
      {
        outputFile << interfaceConc[i][faceI] << tab;
      }
      outputFile << endl;
    }*/

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
