// Required geometrical information
label surfaceID(-1);
forAll (mesh.boundary(), patchI)
{
    if (mesh.boundary()[patchI].name() == "bubble")
    {
        surfaceID = patchI;
    }
}

const polyPatch& bouPatch = mesh.Cf().boundaryField()[surfaceID].patch().patch();

// Get global IDs of faces in boundary patch
labelList patchFaceIDs(bouPatch.size());
forAll (bouPatch,faceI)
{
    patchFaceIDs[faceI] = bouPatch.start()+faceI;
}

// Get IDs of primary adjacent cells
const labelList& adjacentCellIDs = bouPatch.faceCells();

// Get IDs of faces of the cell that are opposite to the patch face
labelList oppFaceIDs(patchFaceIDs.size());

forAll (patchFaceIDs, faceI)
{
    oppFaceIDs[faceI] =
        mesh.cells()[adjacentCellIDs[faceI]].opposingFaceLabel
        (
            patchFaceIDs[faceI],mesh.faces()
        );
}

// Get IDs of secondary adjacent cells
labelList secAdjacentCellIDs(adjacentCellIDs.size());

label globFace = -1;
forAll (oppFaceIDs, faceI)
{
    globFace = oppFaceIDs[faceI];

    if (mesh.owner()[globFace] == adjacentCellIDs[faceI])
    {
        secAdjacentCellIDs[faceI] = mesh.neighbour()[globFace];
    }
    else
    {
        secAdjacentCellIDs[faceI] = mesh.owner()[globFace];
    }
}

// Get IDs of faces of the cell that are opposite to the opponent face
labelList secOppFaceIDs(oppFaceIDs.size());

forAll (oppFaceIDs, faceI)
{
    secOppFaceIDs[faceI] =
        mesh.cells()[secAdjacentCellIDs[faceI]].opposingFaceLabel
        (
            oppFaceIDs[faceI],mesh.faces()
        );
}

// compute distance between boundary face and first cell face opposite
scalarList distF(bouPatch.size(), 0.0);
forAll(bouPatch, faceI)
{
    distF[faceI] = mag(
                   mesh.Cf().boundaryField()[surfaceID][faceI]
                 - mesh.Cf()[oppFaceIDs[faceI]]
               );
}

// area of boundary faces
const scalarField magSf(mag(mesh.Sf().boundaryField()[surfaceID]));
