// create a passive scalar field, T, read BC and read IC
Info << "reading in field T" << endl;
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

Info << "reading in field U" << endl;
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

Info << "reading in transport properties" << endl;

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

Info << "reading in diffusion coefficient" << endl;

dimensionedScalar DT (transportProperties.lookup("DT"));

//#include "createPhi.H" 

Info << "Reading/creating face flux field phi" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf() // [m s-1  * m2]

);
