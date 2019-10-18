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

\*---------------------------------------------------------------------------*/

#include "bubbleSurfaceVelocityComplexFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::
bubbleSurfaceVelocityComplexFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(Zero),
    axis_(Zero),
    normal_(Zero),
    model_name_velocity_(""),
    model_name_radius_("")
{}


Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::
bubbleSurfaceVelocityComplexFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    normal_(dict.lookup("normal")),
    model_name_velocity_(dict.lookupOrDefault<word>("velocity_model", "velocity_model.pt")),
    model_name_radius_(dict.lookupOrDefault<word>("radius_model", "radius_model.pt"))
{
    pyTorch_velocity_ = torch::jit::load(model_name_velocity_);
    pyTorch_radius_ = torch::jit::load(model_name_radius_);
    if (!dict.found("value"))
    {
      updateCoeffs();
    }
}


Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::
bubbleSurfaceVelocityComplexFvPatchVectorField
(
    const bubbleSurfaceVelocityComplexFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    normal_(ptf.normal_),
    model_name_velocity_(ptf.model_name_velocity_),
    model_name_radius_(ptf.model_name_radius_),
    pyTorch_velocity_(ptf.pyTorch_velocity_),
    pyTorch_radius_(ptf.pyTorch_radius_)
{}


Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::
bubbleSurfaceVelocityComplexFvPatchVectorField
(
    const bubbleSurfaceVelocityComplexFvPatchVectorField& rwvpvf
)
:
    fixedValueFvPatchField<vector>(rwvpvf),
    origin_(rwvpvf.origin_),
    axis_(rwvpvf.axis_),
    normal_(rwvpvf.normal_),
    model_name_velocity_(rwvpvf.model_name_velocity_),
    model_name_radius_(rwvpvf.model_name_radius_),
    pyTorch_velocity_(rwvpvf.pyTorch_velocity_),
    pyTorch_radius_(rwvpvf.pyTorch_radius_)
{}


Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::
bubbleSurfaceVelocityComplexFvPatchVectorField
(
    const bubbleSurfaceVelocityComplexFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(rwvpvf, iF),
    origin_(rwvpvf.origin_),
    axis_(rwvpvf.axis_),
    normal_(rwvpvf.normal_),
    model_name_velocity_(rwvpvf.model_name_velocity_),
    model_name_radius_(rwvpvf.model_name_radius_),
    pyTorch_velocity_(rwvpvf.pyTorch_velocity_),
    pyTorch_radius_(rwvpvf.pyTorch_radius_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // compute tangent vector
    const vectorField n(patch().nf());
    const vectorField tau(n ^ normal_);

    // compute polar angle for each patch face
    const vectorField Cf(patch().Cf() - origin_);
    scalar minC = gMin(Cf & axis_);
    const vectorField Cf_min(Cf + mag(minC) * axis_);
    scalarField rf(Cf_min.size(), Zero);

    torch::Tensor phiTensor = torch::ones({Cf.size(), 1}, torch::kFloat64);

    forAll(Cf_min, faceI)
    {
        rf[faceI] = sqrt(Cf_min[faceI] & Cf_min[faceI]);
        phiTensor[faceI][0] = acos((Cf_min[faceI] & axis_) / rf[faceI]);
    }

    // run forward pass to compute tangential velocity
    std::vector<torch::jit::IValue> modelFeatures{phiTensor};
    torch::Tensor uTensor = pyTorch_velocity_.forward(modelFeatures).toTensor();
    torch::Tensor rTensor = pyTorch_radius_.forward(modelFeatures).toTensor();
    auto uAccessor = uTensor.accessor<double,2>();
    auto rAccessor = rTensor.accessor<double,2>();

    vectorField surfaceVelocity(Cf.size(), Zero);
    forAll(surfaceVelocity, faceI)
    {
        scalar dist_inner = mag(rf[faceI] - rAccessor[faceI][0]);
        scalar dist_outer = mag(rf[faceI] - rAccessor[faceI][1]);
        if (dist_inner <= dist_outer)
        {
            surfaceVelocity[faceI] = tau[faceI] * uAccessor[faceI][0];
        }
        else
        {
            surfaceVelocity[faceI] = tau[faceI] * uAccessor[faceI][1];
        }
    }

    vectorField::operator=(surfaceVelocity);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::bubbleSurfaceVelocityComplexFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    os.writeEntry("normal", normal_);
    os.writeEntry("velocity_model", model_name_velocity_);
    os.writeEntry("radius_model", model_name_radius_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        bubbleSurfaceVelocityComplexFvPatchVectorField
    );
}

// ************************************************************************* //
