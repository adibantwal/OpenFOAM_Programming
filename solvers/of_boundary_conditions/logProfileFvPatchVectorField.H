/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::logProfileFvPatchVectorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    for an atmospheric boundary layer profile inlet. This is calcualted from a 
    reference height and reference velocity

Usage
    \table
        Property     | Description                          | Required | Default
        refHeight    | reference height                     | yes |
        refVelocity  | time-based ramping                   | yes  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            logProfile;
        refHeight       100;           
        refVelocity     10;        

    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

    controlDict specification:
    - libs ("liblogProfile.so");

See also
    Foam::fixedValueFvPatchField

SourceFiles
    logProfileFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef logProfileFvPatchVectorField_H
#define logProfileFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class logProfileFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class logProfileFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        scalar refHeight_;
        scalar refVelocity_;

public:

    //- Runtime type information
    TypeName("logProfile");


    // Constructors

        //- Construct from patch and internal field
        logProfileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        logProfileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  logProfileFvPatchVectorField
        //  onto a new patch
        logProfileFvPatchVectorField
        (
            const logProfileFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        logProfileFvPatchVectorField
        (
            const logProfileFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new logProfileFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        logProfileFvPatchVectorField
        (
            const logProfileFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new logProfileFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
