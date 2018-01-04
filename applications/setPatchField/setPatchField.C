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

std::vector<std::vector<std::vector<float> > > getValueArray
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    wordList patchNames,
    unsigned int nTimes,
    std::vector<float>& values
)
{
  word patchName;
  typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
  fieldType field(fieldHeader, mesh);

  unsigned int nPatches = patchNames.size();
  std::vector<std::vector<std::vector<float> > > valueArray;
  valueArray.resize(nPatches);
  for ( unsigned int i = 0 ; i < nTimes ; i++ ) { valueArray[i].resize(nTimes); }
  //TODO: Step through all patches getting their size and setting a data structure
  //forAll(patchNames, patchNameI)
  //{
    //patchName = patchNames[patchNameI];
    //unsigned int patchSize = field.boundaryField()[patchi].size();
  //}
  return valueArray;
}

int setVectorPatchField
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label timei,
    unsigned int nTimes,
    const label patchi,
    std::vector<float>& values,
    int& start_marker,
    bool& done
)
{
    unsigned int nValues = values.size() / 3;

    if (!done)
    {
      typedef GeometricField<vector, fvPatchField, volMesh> fieldType;
      fieldType field(fieldHeader, mesh);

      unsigned int patchSize = field.boundaryField()[patchi].size();
      if ( (start_marker + patchSize) > nValues ) {
        std::cout << "Patch size ( " << patchSize << " ) != number of values ( " << nValues << " ) exiting"<< std::endl;
        return 0;
      }

      Field<vector> patchValues( nValues );

      SubField<vector>
      (
          patchValues,
          patchSize
      ).assign(field.boundaryField()[patchi]);

      int index = 0;
      for (unsigned int i=0; i<patchSize; i++)
      {
          index = start_marker + i*nTimes + timei;
          for(int j = 0; j < 3; ++j) {
            patchValues[i][j] = values[index*3 + j];
          }
      }

      //Info << "Max index " << index << endl;

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
      return nTimes*patchSize;
    }
  return 0;
}

int setScalarPatchField
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label timei,
    unsigned int nTimes,
    const label patchi,
    std::vector<float>& values,
    int& start_marker,
    bool& done
)
{
    if (!done)
    {
      typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;
      fieldType field(fieldHeader, mesh);

      unsigned int patchSize = field.boundaryField()[patchi].size();

      if ( (start_marker + patchSize) > values.size() ) {
        std::cout << "Patch size ( " << patchSize << " ) != number of values ( " << values.size() << " ) exiting"<< std::endl;
        return 0;
      }

      Field<scalar> patchValues( values.size() );

      SubField<scalar>
      (
          patchValues,
          patchSize
      ).assign(field.boundaryField()[patchi]);

      for (unsigned int i=0; i<patchSize; i++)
      {
          int index = start_marker + i*nTimes + timei;
          patchValues[i] = values[index];
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

      return nTimes*patchSize;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    // Read in values
    float f;
    std::vector<float> values;
    while (std::cin >> f) { values.push_back(f); }

    //TODO:: Read in multiple time values and multiple patches in one hit
    // Parse in a huge array from stdin with the values ordered by:
    // cell*times*patch and then parse that into a dictionary or array or arrays
    // which we can then use to set the values
    std::cout << "Read " << values.size() << " values from stdin."<< std::endl;

    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
    argList::addOption("patchNames", "patchNames");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);
    wordList patchNames;

    if (args.optionFound("patchNames")) {
      patchNames = args.optionReadList<word>("patchNames");
    } else {
      patchNames.append(patchName);
    }
    Info << "Calculating for patche(s):" << patchNames << endl;

    unsigned int nTimes = timeDirs.size();
    Info << "Number of times = " << nTimes << endl;

    int counter = 0;
    int start_marker = 0;
    forAll(patchNames, patchNameI)
    {
      start_marker += counter;
      patchName = patchNames[patchNameI];

      Info << "Setting patch " << patchName << endl;
      Info << "Start marker at "<< start_marker << endl;

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
              counter = setScalarPatchField(mesh, io, timeI, nTimes, patchI, values, start_marker, done);
            } else {
              counter = setVectorPatchField(mesh, io, timeI, nTimes, patchI, values, start_marker, done);
            }

            if (!done)
            {
                FatalError
                    << "Only possible to average volVectorFields or volScalarFields."
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
    }
    Info<< "End\n" << endl;

    return 1;
}
