/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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
    fieldManipulation

Description
    This solver shows how to manipulate basic openfoam velocity and pressure fields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

scalar calculatePressure(scalar t, vector x, vector x0, scalar scale);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info <<"Reading in transport properties" << nl << endl;

    IOdictionary transportProperties
    {
        IOobject
        {
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        }
    };
    
    dimensionedScalar nu
    {
        "nu",
        dimViscosity, 
        transportProperties.lookup("nu")
    };

    Info <<"Reading field p\n" <<endl;

    volScalarField p
    {
        IOobject
        {
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        },
        mesh
    };

    Info <<"Reading field U\n" <<endl;

    volVectorField U
    {
        IOobject
        {
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE           
        },
        mesh
    };

    const vector originVector(0.05, 0.05, 0.005);
    const scalar rFarCell = max(
        mag(dimensionedVector("x0",dimLength,originVector)-mesh.C())
    ).value();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        for(label cellI =0; cellI < mesh.C().size(); cellI++)
        {
            p[cellI] = calculatePressure(runTime.time().value(),mesh.C()[cellI], originVector, rFarCell);
        };


        U = fvc::grad(p)*dimensionedScalar("t",dimTime,1.0);
        runTime.write();
    };

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

// Define pressure function

scalar calculatePressure(scalar t, vector x, vector x0, scalar scale)
{
    scalar r (mag(x-x0)/scale);
    scalar rR (1.0/(r + 1e-12));
    scalar f (1.0);

    return Foam::sin(2.0*Foam::constant::mathematical::pi*f*t)*rR;
}
// ************************************************************************* //
