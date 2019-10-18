/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    reactiveSpeciesFoam

Description
    Solves transport equations for passive species with chemical reactions.
    The currently supported reaction types are:
    decayReaction A -kA> P
    (e.g. conversion of species A to P with the reaction rate constant kA)
    singleReaction A+B -kAB> P
    parallelCompetitiveReaction A+B -kAB> P and A+C -kAC> S
    parallelConsecutiveReaction A+B -kAB> P and A+P -kAP> S
    A    - transfer species
    Aphy - physical species transfer (to calculate enhancement)
    B/C  - bulk species
    P    - product
    S    - side product
    The reaction type is defined in the reactionType dictionary.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    if (reaction == "decayReaction")
    {
        #include "solveDecayReaction.H"
    }
    else if(reaction == "singleReaction")
    {
        #include "solveSingleReaction.H"
    }
    else if(reaction == "parallelConsecutiveReaction")
    {
        #include "solveParallelConsecutiveReaction.H"
    }
    else if(reaction == "parallelCompetitiveReaction")
    {
        #include "solveParallelCompetitiveReaction.H"
    }
    else
    {
      Info << "Unknown reaction type " << reaction << endl;
      return 0;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
