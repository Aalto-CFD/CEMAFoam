/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics
    

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "cemaPyjacChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModel(psiReactionThermo);
    makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );



    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        cemaPyjacChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );

}

// ************************************************************************* //
