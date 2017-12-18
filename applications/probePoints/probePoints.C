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

#include "meshSearch.H"
#include "fvCFD.H"
#include "IOprobes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<class FieldType>
void getPatchFaceData
(
    const fvMesh& mesh,
    const IOobject& fieldHeader,
    const label patchI,
    const scalar& zoff,
    List<std::string>& lines,
    bool& done
)
{
    const scalar TOL = 1e-3;
    /*Info << "Trying to read patch " << patchI
      << ", headerClassName:" << fieldHeader.headerClassName()
      << ", FieldTypeName: " << FieldType::typeName << endl;
    */

    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        //Info<< "    Reading " << fieldHeader.headerClassName() << " "
        //    << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);

        Field<typename FieldType::value_type> faceField = field.boundaryField()[patchI].patchInternalField();

        //Set up a stream to write the value to
        std::ostringstream buf;
        OSstream valueString(buf, "value");

	// Now print the actual information
        if (fabs(zoff) < TOL) {
            forAll(faceField, faceI)
            {
                valueString << faceField[faceI];
                if (lines[faceI] != "") { lines[faceI] += ","; }
                lines[faceI] += buf.str();
                buf.str("");buf.clear();
                valueString.flush();
            }
        } else {
            (void)mesh.tetBasePtIs();
            meshSearch meshSearchEngine(mesh);
            const fvPatch& cPatch = mesh.boundary()[patchI];
            const vectorField& faceCenters = cPatch.Cf();
            forAll(faceCenters, faceI)
            {
                point samplePoint = faceCenters[faceI];
                samplePoint[2] += zoff;
                label cellI = meshSearchEngine.findCell(samplePoint);
                if (cellI >= 0) {
                    valueString << field[cellI];
                } else {
                    valueString << faceField[faceI];
                }
                if (lines[faceI] != "") { lines[faceI] += ","; }
                lines[faceI] += buf.str();
                buf.str("");buf.clear();
                valueString.flush();
            }
        }
        done = true;
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noBanner();
    argList::addOption("field", "word", "field to extract patch values");
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

    // Read in the list of probe locations
    IOprobes sniff
    (
        probes::typeName,
        mesh,
        word("probesDict"), // force the use of the system directory
        IOobject::MUST_READ,
        true
    );

    // Get arguments
    word fieldName("");
    args.optionReadIfPresent("field", fieldName);

    label patchI = mesh.boundaryMesh().findPatchID(patchName);
    if (patchI < 0)
    {
        FatalError
            << "Unable to find patch " << patchName << nl
            << exit(FatalError);
    }

    List<std::string> lines(cPatch.size());
    forAll(lines, lineI) { lines[lineI] = ""; }

    if (fieldName != "") {
      forAll(timeDirs, timeI)
      {
        runTime.setTime(timeDirs[timeI], timeI);
        //Info<< "Time = " << runTime.timeName() << endl;

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
            getPatchFaceData<volScalarField>(mesh, io, patchI, zoffset, lines, done);
            getPatchFaceData<volVectorField>(mesh, io, patchI, zoffset, lines, done);

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
    }

   forAll(lines, lineI) { Info << lines[lineI].c_str() << endl; }

    return 0;
}

// ************************************************************************* //
