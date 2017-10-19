/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2011 Mark Pitman
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
    patchRays

Description
    Prints a list of rays to stdout suitable for feeding into the Radiance program
    rtrace in order to get ray-traced irradiance intensity over a patch.
    The output to stdout is in the form:

    xorg yorg zorg xdir ydir zdir (one line per patch face)

    By default rays start originate at the face-cell and point in the direction of the
    face-center.

    The optional argument -out will get rays starting at face-centers and pointing
    outwards in the normal-direction

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noBanner();
    argList::validOptions.insert("awayFromPatch", "");
    argList::validOptions.insert("onlyOrigins", "");
    argList::validArgs.append("patchName");
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

    // Get arguments
    word patchName(args.argRead<word>(1));
    bool outwards = args.optionFound("awayFromPatch");
    bool onlyOrigins = args.optionFound("onlyOrigins");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        //Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        label patchi = mesh.boundaryMesh().findPatchID(patchName);
        if (patchi < 0)
        {
            FatalError
                << "Unable to find patch " << patchName << nl
                << exit(FatalError);
        }

        // Print out Face centers
        const fvPatch& cPatch = mesh.boundary()[patchi];
        const vectorField& faceCenters = cPatch.Cf();

        // Print out cell-center to face-center vector
        Field<Vector<double> > delta;
        int m = 1;
        if (outwards) { delta = cPatch.nf(); m=-1; }
        else { delta = cPatch.delta(); }

        // Now print the actual information
        vector rayOrigin(0.0,0.0,0.0);
        forAll(faceCenters, faceI)
        {
            if (outwards) { rayOrigin = faceCenters[faceI]; }
            else { rayOrigin = (faceCenters[faceI] - delta[faceI]); }

            for (int i=0; i<3; i++) { Info << rayOrigin[i] << " "; }
            if ( not onlyOrigins )
            {
                for (int i=0; i<3; i++) { Info << m*delta[faceI][i] << " "; }
            }

            Info << endl;
        }
    }
    return 0;
}

// ************************************************************************* //
