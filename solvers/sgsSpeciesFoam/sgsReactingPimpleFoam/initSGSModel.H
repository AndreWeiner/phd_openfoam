// Required geometrical information
const fvPatchField<scalar>& pfA = fields[0].boundaryField()[surfaceID];
const polyPatch& bouPatch = pfA.patch().patch();

// Get IDs of faces in boundary patch
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
/*
// Get IDs of third adjacent cells
labelList thirdAdjacentCellIDs(adjacentCellIDs.size());

label globFaceTh = -1;
forAll (secOppFaceIDs, faceI)
{
    globFaceTh = secOppFaceIDs[faceI];

    if (mesh.owner()[globFaceTh] == secAdjacentCellIDs[faceI])
    {
        thirdAdjacentCellIDs[faceI] = mesh.neighbour()[globFaceTh];
    }
    else
    {
        thirdAdjacentCellIDs[faceI] = mesh.owner()[globFaceTh];
    }
}
*/
