/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::bubbleSurfaceVelocitySimpleFvPatchVectorField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition emulating
    the interface velocity at a bubble interface. The velocity is evaluated
    based on the polar angle (axis-symmetry in the azimuth direction assumed)
    and a PyTorch model returning the tangential velocity (scalar product of
    tangent and velocity vector) for a given polar angle.

Usage
    \table
        Property     | Description             | Required    | Default value
        origin       | bubble center of mass in Cartesian coordinates | yes |
        axis         | rise direction          | yes |
        normal       | plane normal vector     | yes |
        model_name   | name of the PyTorch model | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            bubbleSurfaceVelocitySimple;
        origin          (0 0 0);
        axis            (0 1 0);
        normal          (0 0 -1);
        model_name      "bhaga_02_l16_u_tau_ts.pt";
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    bubbleSurfaceVelocitySimpleFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleSurfaceVelocitySimpleFvPatchVectorField_H
#define bubbleSurfaceVelocitySimpleFvPatchVectorField_H

#include <torch/script.h>
#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bubbleSurfaceVelocitySimpleFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bubbleSurfaceVelocitySimpleFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- bubble center of mass
        vector origin_;

        //- rise direction
        vector axis_;

        //- vector normal to 2D mesh / wedge
        vector normal_;

        //- name of the PyTorch model
        word model_name_;

        //- PyTorch model
        torch::jit::script::Module pyTorch_model_;


public:

    //- Runtime type information
    TypeName("bubbleSurfaceVelocitySimple");


    // Constructors

        //- Construct from patch and internal field
        bubbleSurfaceVelocitySimpleFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bubbleSurfaceVelocitySimpleFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given bubbleSurfaceVelocitySimpleFvPatchVectorField
        //  onto a new patch
        bubbleSurfaceVelocitySimpleFvPatchVectorField
        (
            const bubbleSurfaceVelocitySimpleFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        bubbleSurfaceVelocitySimpleFvPatchVectorField
        (
            const bubbleSurfaceVelocitySimpleFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new bubbleSurfaceVelocitySimpleFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bubbleSurfaceVelocitySimpleFvPatchVectorField
        (
            const bubbleSurfaceVelocitySimpleFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
             const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new bubbleSurfaceVelocitySimpleFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
