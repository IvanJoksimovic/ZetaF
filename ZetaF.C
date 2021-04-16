/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "ZetaF.H"
#include "fvOptions.H"

#include "wallFvPatch.H"
#include "wallDist.H"

#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> ZetaF<BasicTurbulenceModel>::Tau() const
{
    return 1.0/omega_;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> ZetaF<BasicTurbulenceModel>::L() const
{
	return CL_*pow(k_,0.5)/omega_;
}

template<class BasicTurbulenceModel>
void ZetaF<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = this->Cmu_*zeta_*k_*Tau();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
ZetaF<BasicTurbulenceModel>::ZetaF
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_		(dimensioned<scalar>::lookupOrAddToDict("Cmu", 			this->coeffDict_, 0.22)),
    COmega2_		(dimensioned<scalar>::lookupOrAddToDict("COmega2", 		this->coeffDict_, 0.9)),
    C1_			(dimensioned<scalar>::lookupOrAddToDict("C1", 			this->coeffDict_, 0.4)),
    C2_			(dimensioned<scalar>::lookupOrAddToDict("C2", 			this->coeffDict_, 0.65)),
    sigmaK_		(dimensioned<scalar>::lookupOrAddToDict("sigmaK", 		this->coeffDict_, 1.1)),
    sigmaOmega_		(dimensioned<scalar>::lookupOrAddToDict("sigmaOmega", 		this->coeffDict_, 1.1)),
    sigmaCDv_		(dimensioned<scalar>::lookupOrAddToDict("sigmaCDv", 		this->coeffDict_, 1.2)),
    sigmaCDt_		(dimensioned<scalar>::lookupOrAddToDict("sigmaCDt", 		this->coeffDict_, 1.6)),
    sigmaZeta_		(dimensioned<scalar>::lookupOrAddToDict("sigmaZeta", 		this->coeffDict_, 1.2)),
    CTau_		(dimensioned<scalar>::lookupOrAddToDict("CTau", 		this->coeffDict_, 6.0)),
    CL_			(dimensioned<scalar>::lookupOrAddToDict("CL", 			this->coeffDict_, 0.36)),
    CEta_		(dimensioned<scalar>::lookupOrAddToDict("CEta", 		this->coeffDict_, 85)),
    Csas_		(dimensioned<scalar>::lookupOrAddToDict("Csas", 		this->coeffDict_, 4)),
    CT2_		(dimensioned<scalar>::lookupOrAddToDict("CT2", 			this->coeffDict_, 1.0)),
    Clim_		(dimensioned<scalar>::lookupOrAddToDict("Clim", 		this->coeffDict_, 0.0)),

    k_ 		(IOobject (IOobject::groupName("k", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    omega_ 	(IOobject (IOobject::groupName("omega", 	alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    epsilon_ 	(IOobject (IOobject::groupName("epsilon", 	alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE ), omega_*k_),
    zeta_ 	(IOobject (IOobject::groupName("zeta", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),

    delta_
    (
        IOobject
        (
            IOobject::groupName("delta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("delta", dimLength, 1e-10)
    ),

    f_ 		(IOobject (IOobject::groupName("f", 		alphaRhoPhi.group()), this->runTime_.timeName(),this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE ), this->mesh_),
    yr_(wallDist::New(this->mesh_).y()),
    mTSmall_("mTSmall", dimensionSet(0, 0, -1, 0, 0, 0, 0),1e-10),
    zetaMin_("zetaMin", dimless, 1e-10),
    fMin_("fMin", dimless/dimTime, 1e-10)

  {
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
    bound(zeta_, zetaMin_);
    bound(f_, fMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
  }


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  template<class BasicTurbulenceModel>
  bool ZetaF<BasicTurbulenceModel>::read()
  {
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        COmega2_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());
        sigmaCDv_.readIfPresent(this->coeffDict());
        sigmaCDt_.readIfPresent(this->coeffDict());
        sigmaZeta_.readIfPresent(this->coeffDict());
        CTau_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        CEta_.readIfPresent(this->coeffDict());
        Csas_.readIfPresent(this->coeffDict());
        CT2_.readIfPresent(this->coeffDict());
        Clim_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
  }


  template<class BasicTurbulenceModel>
  void ZetaF<BasicTurbulenceModel>::correct()
  {
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField G(this->GName(), nut*S2);
    const volScalarField T(this->Tau());
    const volScalarField L(this->L());

    const volScalarField COmega1
    (
        "COmega1",
        (1.4*(1.0 + (0.012/(zeta_+zetaMin_))))-1
    );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAS TERM MIT LIMITER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //dimensionedScalar Psaslim("Psaslim", dimensionSet(0, 0, -2, 0, 0), 0.0);

    //delta_.primitiveFieldRef() = pow(this->mesh().V(), 1.0/3.0);

    //volScalarField LvKdeltaTest("LvKdeltaTest", Clim_*0.11*sqrt(0.41*3.51/(0.8/0.09-0.44))*delta_);
    //volScalarField Lderiv("Lderiv", mag(2.0*symm(fvc::grad(U)))/mag(fvc::laplacian(U)));
    //volScalarField T1lim("T1lim", 40.0*1.775*0.41*mag(2.0*symm(fvc::grad(U)))*1.0/max(Lderiv, LvKdeltaTest)*sqrt(k_));
    //volScalarField ISlim("ISlim", pos(Lderiv-LvKdeltaTest));
    //volScalarField T1_ = 40.0*1.775*0.41*mag(fvc::laplacian(U))*sqrt(k_);

    //volScalarField T2_ = 3.0*k_*max(pow(omega_,-2.0)*(fvc::grad(omega_) & fvc::grad(omega_)), pow(k_,-2.0)*(fvc::grad(k_) & fvc::grad(k_)));
    //volScalarField Psas = max(0.003*(T1_ - CT2_*T2_), Psaslim);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    volScalarField CDv("CDv", (2.0/k_*this->nu()/sigmaCDv_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CDt("CDt", (2.0/k_*nut/sigmaCDt_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CD("CD", CDv+max(CDt,dimensionedScalar("0.0", dimless/dimTime, 0.0)));

    // omega equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
      ==
        alpha*rho*COmega1*G/k_*omega_
      - fvm::SuSp(alpha*rho*COmega2_*omega_, omega_)
      + fvm::Sp(alpha*rho*CD, omega_)
      //+ alpha*rho*Csas_*Psas
    );


    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());

    #include "wallTreatmentOmega.H"

    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::Sp(alpha*rho*omega_, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L), f_)
      - (C1_+(C2_*G/(omega_*k_)))*((zeta_ - 2.0/3.0))/(sqr(L)*T)
    );

    fEqn.ref().relax();
    fvOptions.constrain(fEqn.ref());
    solve(fEqn);
    fvOptions.correct(f_);
    bound(f_, fMin_);

    volScalarField fTilda = f_ - 2.0*this->nu()*zeta_/sqr(yr_);

    // zeta Equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(alpha, rho, zeta_)
      + fvm::div(alphaRhoPhi, zeta_)
      - fvm::laplacian(alpha*rho*DzetaEff(), zeta_)
      ==
        alpha*rho*fTilda
      - fvm::Sp(alpha*rho*G/k_, zeta_)
      + fvOptions(alpha, rho, zeta_)
    );

    zetaEqn.ref().relax();
    fvOptions.constrain(zetaEqn.ref());
    solve(zetaEqn);
    fvOptions.correct(zeta_);
    bound(zeta_, zetaMin_);
    zeta_ = min(zeta_, 2.0);

    epsilon_ = omega_*k_;

    correctNut();

  }



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
