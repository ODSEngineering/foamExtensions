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
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<class FieldType>
void getPatchFaceData
(
    const fvMesh& mesh,
    const meshSearch& meshSearchEngine,
    const IOobject& fieldHeader,
    const label patchI,
    const scalar& zoff,
    List<std::string>& lines,
    bool& done
)
{
    clock_t time_start = std::clock();
    clock_t search_start = std::clock();
    clock_t string_start = std::clock();
    double elapsed_secs = 0;
    double search_time = 0;
    double string_time = 0;

    const scalar TOL = 1e-3;
    /*Info << "Trying to read patch " << patchI
      << ", headerClassName:" << fieldHeader.headerClassName()
      << ", FieldTypeName: " << FieldType::typeName << endl;
    */

    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        Info << "     Reading " << fieldHeader.headerClassName() << " " << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);
        elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
        Info << "     Got FieldType " << elapsed_secs << endl;

        Field<typename FieldType::value_type> faceField = field.boundaryField()[patchI].patchInternalField();
        elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
        Info << "     Got faceField " << elapsed_secs << endl;

        //Set up a stream to write the value to
        std::ostringstream buf;
        OSstream valueString(buf, "value");

        elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
        Info << "     Created string buffer " << elapsed_secs << endl;

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
            elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
            Info << "     Got mesh search engine " << elapsed_secs << endl;
            const fvPatch& cPatch = mesh.boundary()[patchI];
            const vectorField& faceCenters = cPatch.Cf();
            forAll(faceCenters, faceI)
            {
                // Measure the time for the search
                search_start = std::clock();
                point samplePoint = faceCenters[faceI];
                samplePoint[2] += zoff;
                label cellI = meshSearchEngine.findCell(samplePoint);
                // Record the total search time
                search_time += double(std::clock() - search_start) / CLOCKS_PER_SEC;

                // Do the string operations
                string_start = std::clock();
                if (cellI >= 0) {
                    valueString << field[cellI];
                } else {
                    valueString << faceField[faceI];
                }
                if (lines[faceI] != "") { lines[faceI] += ","; }
                lines[faceI] += buf.str();
                buf.str("");buf.clear();
                valueString.flush();
                string_time += double(std::clock() - string_start) / CLOCKS_PER_SEC;
            }
            elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
            Info << "     Probed all sample points ( Search time =" << search_time << ", String time = " << string_time << ") " << elapsed_secs << endl;
        }
        done = true;
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noBanner();
    argList::validArgs.append("patchName");
    argList::addOption("patchNames", "patchNames");
    argList::addOption("field", "word", "field to extract patch values");
    argList::addOption("zOffset", "scalar", "Z-offset above patch");
    argList::addBoolOption("faceData", "include face data (faceCenter, faceNormal & faceArea)");
    argList::addBoolOption("toCSV", "output to CSV file(s)");
#   include "setRootCase.H"

    // Get arguments
    word patchName(args.additionalArgs()[0]);
    word fieldName("");
    scalar zoffset;
    bool inclFaceData = args.optionFound("faceData");
    args.optionReadIfPresent("field", fieldName);
    args.optionReadIfPresent("zOffset", zoffset);
    bool writeToCSV = args.optionFound("toCSV");

    // Get the list of patches to calculate
    wordList patchNames;
    if (args.optionFound("patchNames")) {
      patchNames = args.optionReadList<word>("patchNames");
    } else {
      patchNames.append(patchName);
    }

    if (writeToCSV) {
      Info << "Calculating for patche(s):" << patchNames << endl;
    }

    // Get a cache of patch files
    std::vector<List<std::string>> patchLines(patchNames.count());

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

    // Tetrahedralize and make a MeshSearch object
    (void)mesh.tetBasePtIs();
    meshSearch meshSearchEngine(mesh);

    forAll(patchNames, patchNameI)
    {
      patchName = patchNames[patchNameI];
      if (writeToCSV) {
        Info << "Calculating for patch:" << patchName << endl;
      }

      label patchI = mesh.boundaryMesh().findPatchID(patchName);
      if (patchI < 0)
      {
          FatalError
              << "Unable to find patch " << patchName << nl
              << exit(FatalError);
      }


      const fvPatch& cPatch = mesh.boundary()[patchI];
      const vectorField& faceCenters = cPatch.Cf();
      const vectorField& faceNormals = cPatch.Sf();
      const scalarField& faceAreas = cPatch.magSf();

      //List<std::string> lines(cPatch.size());
      List<std::string> lines = patchLines[patchNameI];
      lines.resize(cPatch.size());
      forAll(lines, lineI) { lines[lineI] = ""; }

      // If requested then get the cell details
      if (inclFaceData) {
          std::ostringstream cellbuf;
          OSstream cellDetails(cellbuf, "celldetails");

          forAll(faceCenters, faceI) {
              cellDetails << faceCenters[faceI] << "," << faceNormals[faceI] << "," << faceAreas[faceI];
              lines[faceI] += cellbuf.str();
              cellbuf.str("");cellbuf.clear();cellDetails.flush();
          }
      }


      forAll(timeDirs, timeI)
      {
        runTime.setTime(timeDirs[timeI], timeI);
        if (writeToCSV) {
          Info << " - Time = " << runTime.timeName() << endl;
        }


        IOobject io
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        if (writeToCSV) { Info << "   - Loaded IOObject complete" << endl; }

        if (io.headerOk())
        {
            mesh.readUpdate();
            if (writeToCSV) { Info << "   - Mesh readUpdate complete" << endl; }

            bool done = false;
            getPatchFaceData<volScalarField>(mesh, meshSearchEngine, io, patchI, zoffset, lines, done);
            getPatchFaceData<volVectorField>(mesh, meshSearchEngine, io, patchI, zoffset, lines, done);
            if (writeToCSV) { Info << "   - getPatchFaceData complete" << endl; }

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

     // Write the results out to a file
     if (writeToCSV) {
       OFstream fWriter ( runTime.path()/(patchName + "." + fieldName) );
       forAll(lines, lineI) { fWriter << lines[lineI].c_str() << endl; }
       Info << "Written patch: " << patchName << " for field: " << fieldName << endl;
     } else {
       forAll(lines, lineI) { Info << lines[lineI].c_str() << endl; }
     }

   }

    return 0;
}

// ************************************************************************* //
