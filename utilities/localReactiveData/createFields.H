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
