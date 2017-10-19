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
    interpPatchField

Description
    Does a vector-interpolation (interpolation in x, y and z components) 
    of the velocity-field between the two nearest angles found.

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "meshSearch.H"
#include "dictionaryEntry.H"
#include "Ostream.H"
#include "fvCFD.H"

vectorField getValues
(
Foam::fvMesh& mesh,
word& fieldName,
label& patchi,
float offset
)
{
    typedef GeometricField<vector, fvPatchField, volMesh> fieldType;
    float TOL = 1e-3;

    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );
    fieldType field0(fieldHeader, mesh);    

    // Get values for the boundary patch
    vectorField U = field0.boundaryField()[patchi].patchInternalField();

    // If there is an offset specified then replace boundary values
    // with values at the offset sampled point
    if (abs(offset) > TOL )
    {
        (void)mesh.tetBasePtIs();
        meshSearch meshSearchEngine(mesh);

        const fvPatch& cPatch = mesh.boundary()[patchi];
        const pointField& faceCenters = cPatch.Cf();

        forAll(faceCenters, faceI)
        {
            point samplePoint = faceCenters[faceI];
            samplePoint[2] += offset;

            label cellI = meshSearchEngine.findCell(samplePoint);
            if (cellI >= 0) {
                U[faceI] = field0[cellI];
            }
        }
    }
    return U;
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");
    argList::validArgs.append("interpTime");
    argList::addOption( "offset", "scalar", "Z-offset above patch" );
#   include "setRootCase.H"

// Avoid printing output by doing this instead:
//#   include "createTime.H"
    Foam::Time runTime
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    instantList timeDirs = timeSelector::select0(runTime, args);

// Also avoid printing out statement of "create mesh for time = X"
//#   include "createMesh.H"
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    word fieldName(args.additionalArgs()[0]);
    word patchName(args.additionalArgs()[1]);
    float interpTime( atof(argv[3]) );
    float offset = 0.0;
    args.optionReadIfPresent("offset", offset);
    
    //Info << "Inputs: " << fieldName << ", " << patchName << ", " << interpTime << ", " << offset << endl;
    //return 1;

    //typedef GeometricField<vector, fvPatchField, volMesh> fieldType;

    // Find the time indices that are either side of the interpTime
    float t0 = -1.0;
    label i0 = 0;
    float t1 = -1.0;
    label i1 = 0;
    float TOL = 1e-3;
    bool onTime = false;

    float minTime = timeDirs[0].value();
    int nTimes = timeDirs.size();
    float maxTime = timeDirs[nTimes-1].value();
    //Info << minTime << ", " << maxTime << ", " << nTimes << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        if ( runTime.value() <= interpTime )
        {
            t0 = runTime.value();
            i0 = timeI;
        }
        else 
        {
            t1 = runTime.value();
            i1 = timeI;
            break;
        }
    }

    if ( i1 <= i0 )
    {
        t0 = maxTime;
        i0 = nTimes-1;
        i1 = 0;
        t1 = minTime;
        t1 += 360.0;
        if (interpTime < t0)
        {
            interpTime += 360.0;
        }
    }

    if ( abs(t0 - interpTime) < TOL )
    {
        onTime = true;
    }

    if ( abs(t1 - interpTime) < TOL )
    {
        i0 = i1;
        t0 = t1;
        onTime = true;
    }

    //Info << "(t0=" << t0 << ",i0=" << i0 << "), (t1=" << t1 << ",i1=" << i1 << "), T = " << interpTime << endl;

    label patchi = mesh.boundaryMesh().findPatchID(patchName);

    // Get the patch values of the first time
    runTime.setTime(timeDirs[i0], i0);
    vectorField U0 = getValues(mesh, fieldName, patchi, offset);

    // Get the patch values of the second time
    vectorField U1 = U0;
    if ( not onTime )
    {
        runTime.setTime(timeDirs[i1], i1);
        U1 = getValues(mesh, fieldName, patchi, offset);
    }

    // Determine the weights for the interpolation
    float p = 0;
    if ( not onTime )
    {
        p = (interpTime - t0)/(t1 - t0);
    }
    
    vector Ui;
    for (int i=0; i<U0.size(); i++)
    {
        Ui = (U1[i] - U0[i])*p + U0[i];
        Info << Ui.component(0) << " " << Ui.component(1) << " " << Ui.component(2) << endl;
    }

    return 1;
}

