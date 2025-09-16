/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::incompressible::RASModels::kEpsH

Description
    Standard k-epsilon turbulence model for incompressible flows.

    The default model coefficients correspond to the following:
    @verbatim
        kEpsHCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            sigmaEps    1.3;
        }
    @endverbatim

SourceFiles
    kEpsH.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsH_H
#define kEpsH_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsH Declaration
\*---------------------------------------------------------------------------*/

class kEpsH
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            //++++++++++++++++++++++++//
            dimensionedScalar Ctau_;
	    dimensionedScalar rSmall_;
            //++++++++++++++++++++++++//
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaEps_;
	    


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_; 
            //+++++++++++++++++++++++// 
	    volScalarField X_;
            volScalarField Y_;
            volScalarField r_;
            volScalarField SinPhi_;
            volScalarField CosPhi_;
	    volTensorField Ltt_;
            volTensorField Rtt_;
            volTensorField LTT_;
            volTensorField RTT_;
            volScalarField Omega_;
            volScalarField TauH_;
            volScalarField OT_;
            volScalarField K1_;
            volScalarField K2_;
            volScalarField K3_;
            volScalarField K4_;
	    volScalarField K5_;
	    volScalarField K6_;
	    volScalarField K7_;
            volTensorField gradU;
            volSymmTensorField nonlinStress_;
            //+++++++++++++++++++++++//

public:

    //- Runtime type information
    TypeName("kEpsH");

    // Constructors

        //- Construct from components
        kEpsH
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~kEpsH()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        
        virtual tmp<volSymmTensorField> R() const;
        
        //+++++++++++++++++++++++++++++++++++++++++++++++++//
        virtual tmp<volTensorField> S() const;
        virtual tmp<volScalarField> SrPhi() const;
        virtual tmp<volScalarField> Srz() const;
        virtual tmp<volScalarField> SPhiz() const;
        virtual tmp<volScalarField> Szz() const;
	virtual tmp<volScalarField> Srr() const;
	virtual tmp<volScalarField> SPhiPhi() const;
        virtual tmp<volTensorField> S_() const;
        virtual tmp<volScalarField> SrPhi_() const;
        virtual tmp<volScalarField> Srz_() const;
        virtual tmp<volScalarField> SPhiz_() const;
        virtual tmp<volScalarField> Szz_() const;
	virtual tmp<volScalarField> Srr_() const;
	virtual tmp<volScalarField> SPhiPhi_() const;
        virtual tmp<volScalarField> Ux_() const;
        virtual tmp<volScalarField> Uy_() const;
        //+++++++++++++++++++++++++++++++++++++++++++++++++//

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
