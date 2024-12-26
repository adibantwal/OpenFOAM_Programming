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
    inputoutput

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H" // This creates a fvMesh system and the instance is called mesh

    const word dictName("randomProperty");

    IOobject dictIO
    {
        dictName, //file name
        mesh.time().constant(), //path to file ... I guess this means the constant folder 
        mesh, // Reference to mesh needed by constructor
        IOobject::MUST_READ

    };

    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn(args.executable()) << "Cannot open custom dictionary "
            << dictName << exit(FatalError);

    dictionary customDict;
    customDict = IOdictionary(dictIO);

    // word userName;
    //customDict.lookup("username") >> userName;

    word userName(customDict.lookupOrDefault<word>("username","NA"));
    scalar offsetValue(customDict.lookupOrDefault<scalar>("offsetValue",5.0));
    bool incompressibleFlag(customDict.lookupOrDefault<bool>("incompressibleFlag",true));
    List<scalar> inputValues(customDict.lookup("inputValues"));
    HashTable<vector,word> inputHashTable(customDict.lookup("inputHashTable"));


    Info << nl << "Looked up the following information: " << nl << nl
        << "Username is: " << userName << nl
        << "Offset value is: " << offsetValue << nl
        << "Incompressible flag is: " << incompressibleFlag << nl
        << "Input values are: " << inputValues << nl
        << "Input hashvalues are: " << inputHashTable << endl;

    // Create output path directoy
    fileName outputDirectory = mesh.time().path()/"postProcessing";
    mkDir(outputDirectory);

    // File pointer to direct the output to
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(outputDirectory/"simulationProperties.dat"));

    // Write to the output file
    outputFilePtr() << "Simulation information" << endl;
    outputFilePtr() << "Simulation run by username: " << userName << endl;

    inputHashTable.insert("finalResult", vector(1,1,1));
    outputFilePtr() << "Final hash table: " << inputHashTable << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
