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
    Foam::bubbleSurfaceVelocityComplexFvPatchVectorField

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
        velocity_model   | name of the PyTorch velocity model | yes |
        radius_model   | name of the PyTorch radius model | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            bubbleSurfaceVelocityComplex;
        origin          (0 0 0);
        axis            (0 1 0);
        normal          (0 0 -1);
        velocity_model  "bhaga_02_l16_u_tau_ts.pt";
        radius_model    "bhaga_02_l16_shape_ts.pt";
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    bubbleSurfaceVelocityComplexFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleSurfaceVelocityComplexFvPatchVectorField_H
#define bubbleSurfaceVelocityComplexFvPatchVectorField_H

#include <torch/script.h>
#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class bubbleSurfaceVelocityComplexFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class bubbleSurfaceVelocityComplexFvPatchVectorField
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

        //- name of the PyTorch velocity model
        word model_name_velocity_;

        //- name of the PyTorch velocity model
        word model_name_radius_;

        //- PyTorch model for the surface velocity
        torch::jit::script::Module pyTorch_velocity_;

        //- PyTorch model for the bubble radius
        torch::jit::script::Module pyTorch_radius_;


public:

    //- Runtime type information
    TypeName("bubbleSurfaceVelocityComplex");


    // Constructors

        //- Construct from patch and internal field
        bubbleSurfaceVelocityComplexFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        bubbleSurfaceVelocityComplexFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given bubbleSurfaceVelocityComplexFvPatchVectorField
        //  onto a new patch
        bubbleSurfaceVelocityComplexFvPatchVectorField
        (
            const bubbleSurfaceVelocityComplexFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        bubbleSurfaceVelocityComplexFvPatchVectorField
        (
            const bubbleSurfaceVelocityComplexFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new bubbleSurfaceVelocityComplexFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        bubbleSurfaceVelocityComplexFvPatchVectorField
        (
            const bubbleSurfaceVelocityComplexFvPatchVectorField&,
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
                new bubbleSurfaceVelocityComplexFvPatchVectorField(*this, iF)
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
