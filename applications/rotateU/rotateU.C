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
    rotateU

Description
    Rotates the U-field velocity vectors according to angle given by time-value (in degrees).

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "quaternion.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    # include "setRootCase.H"

    Foam::Time runTime
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );


    instantList timeDirs = timeSelector::select0(runTime, args);
    bool writeResults = !args.optionFound("noWrite");

    forAll(timeDirs, timeI)
    {
      runTime.setTime(timeDirs[timeI], timeI);

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

      IOobject Uheader
      (
          "U",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ
      );

      if (Uheader.typeHeaderOk<volVectorField>(false))
      {
          Info<< "    Reading U" << endl;
          volVectorField U(Uheader, mesh);

          // Convert to radians
          scalar PI = Foam::constant::mathematical::pi;
          double r = -1.0*runTime.time().value()*PI/180.0;
          quaternion R(vector(0, 0, 1), r);
          tensor T = R.R();
          dimensionedTensor dimT("t", U.dimensions(), T);


          Info<< "    Rotating Velocity Vectors by " << r << " radians" << endl;
          volVectorField Utrans
          (
              IOobject
              (
                  "Utrans",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ
              ),
              U
          );
          transform(Utrans, dimT, Utrans);

          if (writeResults)
          {
              Utrans.write();
          }
      }
      else
      {
          Info<< "    No U" << endl;
      }

      Info<< "\nEnd\n" << endl;
    }
}


// ************************************************************************* //
