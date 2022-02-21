/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::odePyjac

Description
    An ODE solver for chemistry

SourceFiles
    odePyjac.C

\*---------------------------------------------------------------------------*/

#ifndef odePyjac_H
#define odePyjac_H

#include "chemistrySolver.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class odePyjac Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class odePyjac
:
    public chemistrySolver<ChemistryModel>
{
    // Private data

        dictionary coeffsDict_;

        mutable autoPtr<ODESolver> odeSolver_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("odePyjac");


    // Constructors

        //- Construct from thermo
        odePyjac(typename ChemistryModel::reactionThermo& thermo);


    //- Destructor
    virtual ~odePyjac();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "odePyjac.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
