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

\*---------------------------------------------------------------------------*/

#include "odePyjac.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::odePyjac<ChemistryModel>::odePyjac(typename ChemistryModel::reactionThermo& thermo)
:
    chemistrySolver<ChemistryModel>(thermo),
    coeffsDict_(this->subDict("odeCoeffs")),
    odeSolver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::odePyjac<ChemistryModel>::~odePyjac()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ChemistryModel>
void Foam::odePyjac<ChemistryModel>::solve
(
    scalarField& c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cTp_);
    }

    const label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector
    cTp_[0] = p;
    cTp_[1] = T;
    // Update for N-1 species
    for (label i=0; i<nSpecie-1; i++)
    {
        cTp_[i+2] = c[i];
    }

    // We will call ODE solver
    odeSolver_->solve(0, deltaT, cTp_, subDeltaT);

    p = cTp_[0];
    T = cTp_[1];
    scalar csum = 0;

    for (label i=0; i<nSpecie-1; i++)
    {
        c[i] = max(0.0, cTp_[i+2]);
        csum += c[i];
    }
    // Last species
    c[nSpecie-1] = 1.0 - csum;
}


// ************************************************************************* //
