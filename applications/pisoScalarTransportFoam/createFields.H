    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dimensionSet dimSc(0, 0, 0, 0, 0, 0, 0);
    const dimensionSet dimDT(0, 2, -1, 0, 0, 0, 0);

    Info<< "Reading Turbulent Schmidt Number Sc\n";
    /*dimensionedScalar Sc
    (
        transportProperties.lookup("Sc")
    );*/
    dimensionedScalar Sc
    (
        dimensionedScalar::lookupOrDefault
        (
            "Sc",
            transportProperties,
            1,
            dimSc
        )
    );
    Info << Sc << endl;

    Info<< "Reading Molecular Diffusivity DT\n";
    /*dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );*/
    dimensionedScalar DT
    (
        dimensionedScalar::lookupOrDefault
        (
            "DT",
            transportProperties,
            0,
            dimDT
        )
    );
    Info << DT << endl;


#   include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    //IObasicSourceList sources(mesh);
