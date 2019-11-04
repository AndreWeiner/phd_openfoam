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
    speciesFoam

Description
    Solves transport equations for passive species with chemical reactions and
    uses a machine-learning based subgrid-scale model to approximate the
    concentration boundary layer
    Currently supported reaction types are:
    noReaction (only physisorption)
    decayReaction A --k_A-> P
    (e.g. conversion of species A to P with the reaction rate constant kA)
    singleReaction A+B --k_AB-> P
    parallelConsecutiveReaction A+B --k_AB-> P and A+P --k_AP-> S
    A    - transfer species
    B/C  - bulk species
    P    - product
    S    - side product
    The reaction type is defined in the reactionType dictionary.
\*---------------------------------------------------------------------------*/

#include <torch/script.h>
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
    #include "findCellFaceIDs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    if (reaction == "noReaction")
    {
        #include "solvePhysisorption.H"
    }
    else if (reaction == "decayReaction")
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
    else
    {
      Info << "Unknown reaction type " << reaction << endl;
      return 0;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
