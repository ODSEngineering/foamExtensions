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

int main(int argc, char *argv[])
{
    unsigned int i = 0;

    // Read in float-values (one per line)
    float f;
    std::vector<float> values;
    while (std::cin >> f) { values.push_back(f); }

    std::cout << "Read " << values.size() << " values."<< std::endl;

    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    word fieldName(args.argRead<word>(1));
    word patchName(args.argRead<word>(2));

    typedef GeometricField<scalar, fvPatchField, volMesh> fieldType;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field exists
        if (fieldHeader.typeHeaderOk<volScalarField>(false))
        {
            Info<< "    Setting patchField values of "
                << fieldHeader.headerClassName()
                << " " << fieldName << endl;

            fieldType field(fieldHeader, mesh);

            label patchi = mesh.boundaryMesh().findPatchID(patchName);
            if (patchi < 0)
            {
                FatalError
                    << "Unable to find patch " << patchName << nl
                    << exit(FatalError);
            }

            unsigned int patchSize = field.boundaryField()[patchi].size();
            if ( patchSize != values.size() ) { return 0; }

            Field<scalar> patchValues( values.size() );

            SubField<scalar>
            (
                patchValues,
                patchSize
            ) = field.boundaryField()[patchi];

            for (i=0; i<values.size(); i++)
            {
                patchValues[i] = values[i];
            }

            Info<< "    On patch "
                << field.boundaryField()[patchi].patch().name()
                << " set " << patchSize << " values" << endl;

            typename GeometricField<scalar, fvPatchField, volMesh>::
                        Boundary& fieldBf = field.boundaryFieldRef();

            fieldBf[patchi] == SubField<scalar>
            (
                patchValues,
                patchSize
            );

            field.write();

        }
    }

    Info<< "End\n" << endl;

    return 1;
}
