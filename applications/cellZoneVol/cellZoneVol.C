/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    cellZoneVol

Description
    Prints out the total volume of the cellZones in the mesh
    Use as:
    cellZoneVol <cellZoneName>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("cellZoneName");
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
    word cellZoneName(args.additionalArgs()[0]);

    const label cellZoneID = mesh.cellZones().findZoneID(cellZoneName);
    const cellZone& zone = mesh.cellZones()[cellZoneID];
    const cellZoneMesh& zoneMesh = zone.zoneMesh();
    const labelList& cellsZone = zoneMesh[cellZoneID]; //list of all radiatorZone cell ID's

    scalar cellZoneVol(0);
    forAll(cellsZone, cI)
    {
    cellZoneVol += mesh.V()[cellsZone[cI]];
    }

    reduce(cellZoneVol,sumOp<scalar>());
    Info << cellZoneVol << nl << endl;
    return 0;
}



// ************************************************************************* //
