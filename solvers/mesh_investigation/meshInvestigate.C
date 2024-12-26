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
    meshhInvestigate

Description
    This solver investigates the mesh
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "The most recent time folder is " << runTime.timeName() << nl
    << "The mesh has " << mesh.C().size() << " cells and " <<mesh.Cf().size()
    << " internal faces." <<endl;

    Info << endl;

    // Accessing the internal cells
    for (label faceI = 0; faceI < mesh.owner().size() ; faceI++)
        if (faceI%100 == 0)
            Info << "Internal face " << faceI << " with center at " << mesh.Cf()[faceI]
             << " belongs to cell " << mesh.owner()[faceI]
             << " and has a neighbor " << mesh.neighbour()[faceI] << endl;

    Info << endl;

    // Accessing the boundaries
    forAll(mesh.boundaryMesh(),patchI)
        Info << "Patch " << patchI << " = " << mesh.boundaryMesh()[patchI].name() << " with "
        << mesh.boundary()[patchI].Sf().size() << " faces. This starts at face "
        << mesh.boundary()[patchI].start()
        << ". The surface area of the patch " <<  mesh.boundaryMesh()[patchI].name() << " is "
        << gSum(mesh.boundary()[patchI].magSf()) << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
