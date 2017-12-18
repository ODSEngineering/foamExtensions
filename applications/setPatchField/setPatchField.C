/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    setPatchField

Description
    Reads scalar values line-by-line from stdin and sets a non-uniform value
    in a SCALAR PATCH FIELD.

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "dictionaryEntry.H"
#include "Ostream.H"
#include "fvCFD.H"

std::vector<float> lineToVect(std::string str, char delimiter) {
  std::vector<float> vect;
  std::stringstream ss(str); // Turn the string into a stream.
  string tok;

  while(getline(ss, tok, delimiter)) {
    vect.push_back( atof(tok.c_str()) );
  }

  return vect;
}


void setVectorPatchField
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label patchi,
    std::vector<float>& values,
    bool& done
)
{
    unsigned int nValues = values.size() / 3;

    if (!done)
    {
      typedef GeometricField<vector, fvPatchField, volMesh> fieldType;
      fieldType field(fieldHeader, mesh);

      unsigned int patchSize = field.boundaryField()[patchi].size();
      if ( patchSize != nValues ) {
        std::cout << "Patch size ( " << patchSize << " ) != number of values ( " << nValues << " ) exiting"<< std::endl;
        return;
      }

      Field<vector> patchValues( nValues );

      SubField<vector>
      (
          patchValues,
          patchSize
      ).assign(field.boundaryField()[patchi]);

      for (unsigned int i=0; i<nValues; i++)
      {
          for(int j = 0; j < 3; ++j) {
            patchValues[i][j] = values[i*3 + j];
          }
      }

      Info<< "    On patch "
          << field.boundaryField()[patchi].patch().name()
          << " set " << patchSize << " values" << endl;

      field.boundaryField()[patchi] == SubField<vector>
      (
          patchValues,
          patchSize
      );

      field.write();

      done = true;
    }
}

void setScalarPatchField
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label patchi,
    std::vector<float>& values,
    bool& done
)
{
    if (!done)
    {
      typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
      fieldType field(fieldHeader, mesh);

      unsigned int patchSize = field.boundaryField()[patchi].size();
      if ( patchSize != values.size() ) {
        std::cout << "Patch size ( " << patchSize << " ) != number of values ( " << values.size() << " ) exiting"<< std::endl;
        return;
      }

      Field<scalar> patchValues( values.size() );

      SubField<scalar>
      (
          patchValues,
          patchSize
      ).assign(field.boundaryField()[patchi]);

      for (unsigned int i=0; i<values.size(); i++)
      {
          patchValues[i] = values[i];
      }

      Info<< "    On patch "
          << field.boundaryField()[patchi].patch().name()
          << " set " << patchSize << " values" << endl;

      field.boundaryField()[patchi] == SubField<scalar>
      (
          patchValues,
          patchSize
      );

      field.write();

      done = true;
    }
}

int main(int argc, char *argv[])
{
    // Read in values
    float f;
    std::vector<float> values;
    while (std::cin >> f) { values.push_back(f); }

    std::cout << "Read " << values.size() << " values from stdin."<< std::endl;

    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);

    // Check if this patch exists
    label patchI = mesh.boundaryMesh().findPatchID(patchName);
    if (patchI < 0)
    {
        FatalError
            << "Unable to find patch " << patchName << nl
            << exit(FatalError);
    }

    forAll(timeDirs, timeI)
    {
      runTime.setTime(timeDirs[timeI], timeI);
      Info<< "Time = " << runTime.timeName() << endl;

      IOobject io
      (
          fieldName,
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
      );


      if (io.headerOk())
      {
          mesh.readUpdate();

          bool done = false;
          if ( io.headerClassName() == "volScalarField") {
            setScalarPatchField(mesh, io, patchI, values, done);
          } else {
            setVectorPatchField(mesh, io, patchI, values, done);
          }

          if (!done)
          {
              FatalError
                  << "Only possible to average volFields."
                  << " Field " << fieldName << " is of type "
                  << io.headerClassName()
                  << nl << exit(FatalError);
          }
      }
      else
      {
          Info<< "    No field " << fieldName << endl;
      }
    }

    Info<< "End\n" << endl;

    return 1;
}
