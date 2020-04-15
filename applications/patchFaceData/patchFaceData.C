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
    const wordList& patchNames,
    const scalar& zoff,
    List<List<std::string> >& patchLines,
    List<List<point> >& patchSamplePoints,
    List<List<label> >& patchSampleCells,
    bool& done,
    bool& writeToCSV
)
{
    //clock_t time_start = std::clock();
    //clock_t search_start = std::clock();
    //clock_t string_start = std::clock();
    //double elapsed_secs = 0;
    //double search_time = 0;
    //double string_time = 0;

    const scalar TOL = 1e-3;

    if (!done && fieldHeader.headerClassName() == FieldType::typeName)
    {
        //Info << "     Reading " << fieldHeader.headerClassName() << " " << fieldHeader.name() << endl;

        FieldType field(fieldHeader, mesh);
        //elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
        //Info << "     Got FieldType " << elapsed_secs << endl;

        forAll(patchNames, patchNameI)
        {
          word patchName = patchNames[patchNameI];
          label patchI = mesh.boundaryMesh().findPatchID(patchName);
          if (writeToCSV) {
            Info << "     Getting data from patch " <<  patchName << " (size " << patchLines[patchNameI].size() << ")" << endl;
          }


          Field<typename FieldType::value_type> faceField = field.boundaryField()[patchI].patchInternalField();
          //elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
          //Info << "     Got faceField " << elapsed_secs << endl;

          //Set up a stream to write the value to
          std::ostringstream buf;
          OSstream valueString(buf, "value");

          //elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
          //Info << "     Created string buffer " << elapsed_secs << endl;

  	      // Now print the actual information
          if (fabs(zoff) < TOL) {
              forAll(faceField, faceI)
              {
                  valueString << faceField[faceI];
                  if (patchLines[patchNameI][faceI] != "") { patchLines[patchNameI][faceI] += ","; }
                  patchLines[patchNameI][faceI] += buf.str();
                  buf.str("");buf.clear();
                  valueString.flush();
              }
          } else {
              //elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
              //Info << "     Got mesh search engine " << elapsed_secs << endl;
              const fvPatch& cPatch = mesh.boundary()[patchI];
              const vectorField& faceCenters = cPatch.Cf();
              forAll(faceCenters, faceI)
              {
                  //point samplePoint = patchSamplePoints[patchNameI][faceI];
                  label cellI = patchSampleCells[patchNameI][faceI];

                  // Do the string operations
                  //string_start = std::clock();
                  if (cellI >= 0) {
                      valueString << field[cellI];
                  } else {
                      // Return the surface value if the sample point is not found
                      valueString << faceField[faceI];
                  }
                  if (patchLines[patchNameI][faceI] != "") { patchLines[patchNameI][faceI] += ","; }
                  patchLines[patchNameI][faceI] += buf.str();
                  buf.str("");buf.clear();
                  valueString.flush();
                  //string_time += double(std::clock() - string_start) / CLOCKS_PER_SEC;
              }
              //elapsed_secs = double(std::clock() - time_start) / CLOCKS_PER_SEC;
              //Info << "     Probed all sample points ( Search time =" << search_time << ", String time = " << string_time << ") " << elapsed_secs << endl;
          }
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
    List<List<std::string> > patchLines(patchNames.size());

    // Get the points that we will sample on each patch
    List<List<point> > patchSamplePoints(patchNames.size());
    List<List<label> > patchSampleCells(patchNames.size());

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

    // First go through each patch and get the patch details
    forAll(patchNames, patchNameI)
    {
      patchName = patchNames[patchNameI];
      if (writeToCSV) {
        Info << "Getting face details for patch: " << patchName << endl;
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
      patchLines[patchNameI].resize(cPatch.size());
      forAll(patchLines[patchNameI], lineI) { patchLines[patchNameI][lineI] = ""; }

      // If requested then get the cell details
      if (inclFaceData) {
          std::ostringstream cellbuf;
          OSstream cellDetails(cellbuf, "celldetails");

          forAll(faceCenters, faceI) {
              cellDetails << faceCenters[faceI] << "," << faceNormals[faceI] << "," << faceAreas[faceI];
              patchLines[patchNameI][faceI] += cellbuf.str();
              cellbuf.str("");cellbuf.clear();cellDetails.flush();
          }
      }

      // Assume a static mesh so get the sample points out now
      patchSamplePoints[patchNameI].resize(cPatch.size());
      patchSampleCells[patchNameI].resize(cPatch.size());
      // Get the values for all 
      forAll(faceCenters, faceI) {
        point samplePoint = faceCenters[faceI];
        samplePoint[2] += zoffset;
        label cellI = meshSearchEngine.findCell(samplePoint);

        patchSamplePoints[patchNameI][faceI] = samplePoint;
        patchSampleCells[patchNameI][faceI] = cellI;
      }
    }



    // Now step through each time and extract the probe data
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
      //if (writeToCSV) { Info << "   - Loaded IOObject complete" << endl; }

      if (io.headerOk())
      {
          mesh.readUpdate();
          //if (writeToCSV) { Info << "   - Mesh readUpdate complete" << endl; }

          // For each patch extract the patch probe data
          bool done = false;
          getPatchFaceData<volScalarField>(mesh, meshSearchEngine, io, patchNames, zoffset, patchLines, patchSamplePoints, patchSampleCells, done, writeToCSV);
          getPatchFaceData<volVectorField>(mesh, meshSearchEngine, io, patchNames, zoffset, patchLines, patchSamplePoints, patchSampleCells, done, writeToCSV);

          //if (writeToCSV) { Info << "   - getPatchFaceData complete" << endl; }

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
   forAll(patchNames, patchNameI)
   {
     if (writeToCSV) {
       OFstream fWriter ( runTime.path()/(patchNames[patchNameI] + "." + fieldName) );
       forAll(patchLines[patchNameI], lineI) { fWriter << patchLines[patchNameI][lineI].c_str() << endl; }
       Info << "Written patch: " << patchNames[patchNameI] << " for field: " << fieldName << endl;
     } else {
       forAll(patchLines[patchNameI], lineI) { Info << patchLines[patchNameI][lineI].c_str() << endl; }
     }
   }

    return 0;
}

// ************************************************************************* //
