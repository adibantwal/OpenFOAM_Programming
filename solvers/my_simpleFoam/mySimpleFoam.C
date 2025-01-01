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
    mySimpleFoam

Description
    Implement a light-weight SIMPLE algorithm in OpenFOAM. Optimization and 
    convergence checks are neglected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Relaxation factor
    scalar alpha;
    fvSolution.lookup("alpha") >> alpha;
    scalar pRefCell;
    fvSolution.lookup("pRefCell") >> pRefCell;
    scalar pRefValue;
    fvSolution.lookup("pRefValue") >> pRefValue;

    Info << "Read the following parameters: " << endl;
    Info << "Relaxation factor: " << alpha << endl; 
    Info << "Index of cell for reference pressure: " << pRefCell << endl; 
    Info << "Value of reference pressure: " << pRefValue << endl; 

    while (runTime.loop())
    {
        Info << nl << "Iteration: " << runTime.timeName() << endl;

        //Define the momentum eqn.
        fvVectorMatrix UEqn
        (
            fvm::div(phi,U) - fvm::laplacian(nu,U) == -fvc::grad(p)

        );

        //Solve momentum equation for current value of pressure
        UEqn.solve();

        volScalarField A = UEqn.A();
        volVectorField H = UEqn.H();

        volScalarField A_inv = 1.0/A;
        surfaceScalarField A_inv_flux = fvc::interpolate(A_inv);
        volVectorField HbyA = A_inv*H;

        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv_flux,p) == fvc::div(HbyA)
        );

        //Set reference pressure for the equation
        pEqn.setReference(pRefCell,pRefValue);

        pEqn.solve(); 

        //Explicit under-relaxation of pressure equation
        p = alpha*p + (1-alpha)*p_old;

        U = (A_inv*H) - (A_inv*fvc::grad(p));

        phi = fvc::interpolate(U) & mesh.Sf();

        //Update boundary conditions
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        p_old = p;

        runTime.write();

    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
