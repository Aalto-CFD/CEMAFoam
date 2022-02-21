/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "cemaPyjacChemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::cemaPyjacChemistryModel
(
    ReactionThermo& thermo
)
:
    BasicChemistryModel<ReactionThermo>(thermo),
    ODESystem(),
    Y_(this->thermo().composition().Y()),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_
    (
        BasicChemistryModel<ReactionThermo>::template getOrDefault<scalar>
        (
            "Treact",
            0.0
        )
    ),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    sp_enthalpy_(nSpecie_),
    nElements_(BasicChemistryModel<ReactionThermo>::template get<label>("nElements")),
    chemJacobian_(nSpecie_),
    cem_
    (
        IOobject
        (
            "cem",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("cem", dimless, 0),
        extrapolatedCalculatedFvPatchScalarField::typeName
    )
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
            )
        );
    }

    Info<< "cemaPyjacChemistryModel: Number of species = " << nSpecie_
        << "\n  and reactions (from reaction file, expected 0 with PyJac) = " << nReaction_ << endl; 
        // Note that nReaction_ should be updated with PyJAC
        // PERHAPS TO OVERWRITE IN THE SRC DURING DYNAMIC BINDING
    Info<< "cemaPyjacChemistryModel: Number of elements = " << nElements_ << endl; 
 
    if (this->chemistry_) {
        Info << "\n Evaluating species enthalpy of formation using PyJac\n" << endl;
        //- Enthalpy of formation for all species
        std::vector<scalar> sp_enth_form(nSpecie_, 0.0);
        //- Enthalpy of formation is taken from pyJac at T-standard (chem_utils.h)
        eval_h(298.15, sp_enth_form.data());
        for (label i = 0; i < nSpecie_; ++i)
        { 
            sp_enthalpy_[i] = sp_enth_form[i];
        }
        chemJacobian_ = Zero;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::
~cemaPyjacChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    dcdt = Zero;

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            const label si = R.lhs()[s].index;
            const scalar sl = R.lhs()[s].stoichCoeff;
            dcdt[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            const label si = R.rhs()[s].index;
            const scalar sr = R.rhs()[s].stoichCoeff;
            dcdt[si] += sr*omegai;
        }
    }
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const Reaction<ThermoType>& R = reactions_[index];
    scalar w = omega(R, c, T, p, pf, cf, lRef, pr, cr, rRef);
    return(w);
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c);
    const scalar kr = R.kr(kf, p, T, c);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0.0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0.0), exp);
        }
    }
    cf = max(c[lRef], 0.0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0.0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0.0), exp);
        }
    }
    cr = max(c[rRef], 0.0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo, class ThermoType>
void Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    std::vector<double> TY(nSpecie_+1, 0.0);
    // if TY has N+1 elements, diff(TY) has N elements
    std::vector<double> dTYdt(nSpecie_, 0.0);

    const scalar p = c[0];
    const scalar T = c[1];

    scalar csum = 0.0;
    forAll(c_, i)
    {
        c_[i] = max(c[i+2], 0.0);
        csum += c_[i];
    }
    // Then we exclude last species from csum and dump all residuals
    // into last species to ensure mass conservation
    csum -= c_[nSpecie_-1];
    c_[nSpecie_-1] = 1.0 - csum;

    TY[0] = T;
    forAll(c_, i)
    { 
        TY[i+1] = c_[i];
    }

    dydt(0, p, TY.data(), dTYdt.data());

    // dp/dt = 0
    dcdt[0] = 0.0;

    // Back substitute into dcdt (dcdt has nSpecie+1 elements for diff(PTY))
    for (label i = 0; i < nSpecie_; ++i)
    {
        dcdt[i+1] = dTYdt[i];
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    std::vector<double> TY(nSpecie_+1, 0.0); //###
    std::vector<double> dfdy(nSpecie_*nSpecie_, 0.0); //###

    const scalar p = c[0];
    const scalar T = c[1];

    scalar csum = 0.0;
    forAll(c_, i)
    {
        c_[i] = max(c[i+2], 0.0);
        csum += c_[i];
    }
    // Then we exclude last species from csum and instead dump all
    // residuals into last species to ensure mass conservation
    csum -= c_[nSpecie_-1];
    c_[nSpecie_-1] = 1.0 - csum;

    dfdc = Zero;

    TY[0] = T;
    // Assign nSpecies-1 species mass fractions to the TY vector
    forAll(c_, i)
    {
        TY[i+1] = c_[i];
    }

    eval_jacob(0, p, TY.data(), dfdy.data());

    // Back substitution to update dfdc

    // Assign first row and column to zero as they correspond to const pressure
    for (label j = 0; j < nSpecie_ + 1; ++j)
    {
        dfdc(0,j) = 0.0;
        dfdc(j,0) = 0.0;
    }

    label k = 0;
    // Loop cols
    for (label j = 1; j < nSpecie_+1; ++j)
    {
        // Loop rows
        for (label i = 1; i < nSpecie_+1; ++i)
        {
            dfdc(i,j) = dfdy[k + i - 1];
            chemJacobian_(i-1,j-1) = dfdy[k + i - 1];
        }
        k += nSpecie_;
    }

    // Note that dcdt is not needed in most ODE solvers so here we just return 0
    dcdt = Zero;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("small", dimTime, SMALL),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    const label nReaction = reactions_.size();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];

            scalar cSum = 0.0;

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                cSum += c_[i];
            }

            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];

                omega(R, c_, Ti, pi, pf, cf, lRef, pr, cr, rRef);

                forAll(R.rhs(), s)
                {
                    tc[celli] += R.rhs()[s].stoichCoeff*pf*cf;
                }
            }

            tc[celli] = nReaction*cSum/tc[celli];
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy / dimVolume / dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                // const scalar hi = specieThermo_[i].Hc();
                scalar hi = sp_enthalpy_[i];
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermo_[i].W();
        }

        const scalar w = omegaI
        (
            ri,
            c_,
            Ti,
            pi,
            pf,
            cf,
            lRef,
            pr,
            cr,
            rRef
        );

        RR[celli] = w*specieThermo_[si].W();
    }

    return tRR;
}


template<class ReactionThermo, class ThermoType>
void Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermo_[i].W();
        }

        omega(c_, Ti, pi, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt_[i]*specieThermo_[i].W();
        }
    }
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                // c_[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                c_[i] = Y_[i][celli];
                c0[i] = c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                // Calls ode::solve() from chemistrySolver
                this->solve(c_, Ti, pi, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<nSpecie_; i++)
            {
                // CHEMICAL SOURCE TERM PER SPECIES
                // (c_[i] - c0[i])*specieThermo_[i].W()/deltaT[celli]; // ###
                this->RR_[i][celli] = rhoi*(this->c_[i] - c0[i])/deltaT[celli];
            }

            if (this->time().write()) {
                    // Info << "\nCELL: " << celli << "\t Temperature = " << Ti << endl;
                    scalar cem_cell;
                    cema(cem_cell);
                    cem_[celli] = cem_cell;
            }
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }

    }

    cem_.correctBoundaryConditions();

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}

template<class ReactionThermo, class ThermoType>
void Foam::cemaPyjacChemistryModel<ReactionThermo, ThermoType>::cema
(
    scalar& cem
) const
{

    const Foam::EigenMatrix<scalar> EM(chemJacobian_);
    DiagonalMatrix<scalar> EValsRe(EM.EValsRe());
    DiagonalMatrix<scalar> EValsIm(EM.EValsIm());

    DiagonalMatrix<scalar> EValsMag(EValsRe.size(), 0.0);
    forAll(EValsRe, i)
    {
        EValsMag[i] = (EValsRe[i]*EValsRe[i] + EValsIm[i]*EValsIm[i]);
    }

    // Sort eigenvalues in ascending order, and track indices
    const auto ascend = [&](scalar a, scalar b){ return a < b; };
    const List<label> permut(EValsMag.sortPermutation(ascend));

    // Skip conservation modes for elements and temperature
    for (label i=0; i<nElements_+1; ++i) 
    {
        label idx = permut[i];
        EValsRe[idx] = -1E30;
    }

    cem = gMax(EValsRe);
}

// ************************************************************************* //
