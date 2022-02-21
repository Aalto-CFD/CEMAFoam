/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation  -- Author: Mahmoud Gadalla
-------------------------------------------------------------------------------
License
    This file is part of CEMAFoam, derived from OpenFOAM.

    https://github.com/Aalto-CFD/CEMAFoam

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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
inline Foam::label
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::nEqns() const
{
    // nEqns = number of species (N-1) + temperature + pressure
    return nSpecie_ + 1;
}


template<class ReactionThermo, class ThermoType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::RR()
{
    return RR_;
}


template<class ReactionThermo, class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType>>&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::reactions() const
{
    return reactions_;
}


template<class ReactionThermo, class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::specieThermo() const
{
    return specieThermo_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::label
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::label
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalar
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::Treact() const
{
    return Treact_;
}


template<class ReactionThermo, class ThermoType>
inline Foam::scalar&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::Treact()
{
    return Treact_;
}


template<class ReactionThermo, class ThermoType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::RR
(
    const label i
) const
{
    return RR_[i];
}

template<class ReactionThermo, class ThermoType>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::RR
(
    const label i
)
{
    return RR_[i];
}


// ************************************************************************* //
