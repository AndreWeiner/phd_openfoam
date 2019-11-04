// phi is necessary to correct some of the boundary conditions
surfaceScalarField phi
(
    IOobject
    (
      "phi",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
);

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
   fields[i].correctBoundaryConditions();
}
