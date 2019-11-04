// determine patch ID
label surfaceID(-1);
forAll (mesh.boundary(), patchI)
{
    if (mesh.boundary()[patchI].name() == "bubble")
    {
        surfaceID = patchI;
    }
}

// compute surface unit normal vector
vectorField surfaceNormal = - mesh.Sf().boundaryField()[surfaceID]/
                              mag(mesh.Sf().boundaryField()[surfaceID]);

// face centers of interface patch
vectorField surfaceCf = mesh.Cf().boundaryField()[surfaceID];

// get IDs of faces in boundary patch
const polyPatch& interfacePatch = mesh.Sf().boundaryField()[surfaceID].patch().patch();
labelList boundaryCellIDs = interfacePatch.faceCells();
std::vector<std::vector<label>> faceIDs(boundaryCellIDs.size());
std::vector<std::vector<label>> cellIDs(boundaryCellIDs.size());
std::vector<std::vector<scalar>> distCf(boundaryCellIDs.size());

forAll(interfacePatch, faceI)
{
    std::vector<label> faces(1, interfacePatch.start() + faceI);
    std::vector<label> cells(1, boundaryCellIDs[faceI]);
    faceIDs[faceI] = faces;
    cellIDs[faceI] = cells;
}

forAll(interfacePatch, faceI)
{
    std::vector<label> faces = faceIDs[faceI];
    std::vector<label> cells = cellIDs[faceI];

    std::vector<scalar> curDist(1, 0.0);
    do {
        // find ID of opposite face
        label lastFaceID = faces.back();
        label lastCellID = cells.back();
        label nextFaceID = mesh.cells()[lastCellID].opposingFaceLabel(
            lastFaceID, mesh.faces()
        );
        // find ID of next cell
        if (mesh.owner()[nextFaceID] == lastCellID)
        {
            cells.push_back(mesh.neighbour()[nextFaceID]);
        }
        else
        {
            cells.push_back(mesh.owner()[nextFaceID]);
        }
        faces.push_back(nextFaceID);
        // compute distance from boundary
        curDist.push_back(mag(mesh.Cf()[nextFaceID] - mesh.Cf()[faces.front()]));
    } while (curDist.back() < maxDist);
    faceIDs[faceI] = faces;
    cellIDs[faceI] = cells;
    distCf[faceI] = curDist;
}
