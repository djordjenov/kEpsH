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

\*---------------------------------------------------------------------------*/

#include "kEpsH.H"
#include "addToRunTimeSelectionTable.H"


#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsH, 0);
addToRunTimeSelectionTable(RASModel, kEpsH, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsH::kEpsH
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//

    Ctau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau",
            coeffDict_,
            0.05
        )
    ),
     
    rSmall_
    (
        dimensionedScalar
        (
            "rSmall",
            dimensionSet (0,1,0,0,0,0,0),
            0.000001
        )
    ),

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),

    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),

    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//
    X_                                               //declared in kEpsH.H
    ( 
          mesh_.C().component(vector::X)             
    ),

    Y_                                               
    (
          mesh_.C().component(vector::Y)
    ),

    r_                                               
    (
          sqrt(sqr(X_)+sqr(Y_))+rSmall_
    ),

    SinPhi_                                        
    (
          Y_/r_
    ),

    CosPhi_                                        
    (
          X_/r_
    ),

    Ltt_  //Left transformation tensor 
    (
        CosPhi_*OneOne(tensor::one)+SinPhi_*OneTwo(tensor::one)
       -SinPhi_*TwoOne(tensor::one)+CosPhi_*TwoTwo(tensor::one)
       + ThreeThree(tensor::one) 
    ),

    Rtt_  //Right transformation tensor 
    (
        CosPhi_*OneOne(tensor::one)-SinPhi_*OneTwo(tensor::one)
       +SinPhi_*TwoOne(tensor::one)+CosPhi_*TwoTwo(tensor::one) 
       + ThreeThree(tensor::one) 
    ),

    LTT_  //Left transformation tensor 
    (
        CosPhi_*OneOne(tensor::one)-SinPhi_*OneTwo(tensor::one)
       +SinPhi_*TwoOne(tensor::one)+CosPhi_*TwoTwo(tensor::one) 
       + ThreeThree(tensor::one) 
    ),

    RTT_  //Right transformation tensor 
    (
        CosPhi_*OneOne(tensor::one)+SinPhi_*OneTwo(tensor::one)
       -SinPhi_*TwoOne(tensor::one)+CosPhi_*TwoTwo(tensor::one) 
       + ThreeThree(tensor::one) 
    ),

    Omega_
    (
            IOobject
            (
                "Omega",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mag(-Ux_()*SinPhi_+Uy_()*CosPhi_)/r_      //Ugaona brzina vihora
    ),

    TauH_
    (
            IOobject
            (
                "TauH",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Ctau_*k_/(epsilon_ + epsilonSmall_)   //Hambin vremenski razmer
    ),

    OT_
     (
            IOobject
            (
                "OT",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Omega_*TauH_                          //Pomocni kontrolni parametar
     ),

    K1_                                             
    (
            IOobject
            (
                "K1",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            4*(OT_)/(1+4*sqr(OT_))
    ),

    K2_                                           
    (
            IOobject
            (
                "K2",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
           2/(1+4*sqr(OT_))
    ),

    K3_                                             
    (
            IOobject
            (
                "K3",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            2/(1+sqr(OT_))   
    ),

    K4_                                            
    (
            IOobject
            (
                "K4",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            2*(OT_)/(1+sqr(OT_))
    ),

    K5_
    (
            IOobject
            (
                "K5",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            2*(OT_)/(1+4*sqr(OT_))
    ),

    K6_
    (
            IOobject
            (
                "K6",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            4*sqr(OT_)/(1+4*sqr(OT_))
    ),

    K7_
    (
            IOobject
            (
                "K7",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            2*(1+2*sqr(OT_))/(1+4*sqr(OT_))
    ),


    gradU(fvc::grad(U_)), // DjN, 22/Aug/2018

    nonlinStress_
    (
     "nonlinStress",
   symm(nut_*symm( //1.L 
		   LTT_&
		   ( //2.L
		   ( //3.L
		   ((
		   -K1_*OneOne(tensor::one)
		   -K2_*OneTwo(tensor::one)
		   -K2_*TwoOne(tensor::one)
		   +K1_*TwoTwo(tensor::one)
		   )
		   &
		   (
		   SrPhi_()*I
		   ))
		   +
		   ((
		   -K3_*OneThree(tensor::one)
		   +K4_*TwoThree(tensor::one)
		   -K3_*ThreeOne(tensor::one)
		   +K4_*ThreeTwo(tensor::one)
		   )
	 	   &
		   (
		   Srz_()*I
		   ))
		   +
		   ((
		   -K4_*OneThree(tensor::one)
		   -K3_*TwoThree(tensor::one)
		   -K4_*ThreeOne(tensor::one)
		   -K3_*ThreeTwo(tensor::one)
		   )
		   &
		   (
		   SPhiz_()*I
		   ))
		   +
		   ((
		   -2*ThreeThree(tensor::one)
		   )
		   &
		   (
		   Szz_()*ThreeThree(tensor::one)
		   ))
		   +
                   ((
                   -K7_*OneOne(tensor::one)
                   +K5_*OneTwo(tensor::one)
                   +K5_*TwoOne(tensor::one)
                   -K6_*TwoTwo(tensor::one)
                   ) 
                   &
                   ( 
                   Srr_()*I
                   ))
                   +
                   ((
                   -K6_*OneOne(tensor::one)
                   -K5_*OneTwo(tensor::one)
                   -K5_*TwoOne(tensor::one)
                   -K7_*TwoTwo(tensor::one)
                   )
                   &
                   (
                   SPhiPhi_()*I
                   ))
	 	   ) //3.R
		   ) //2.R
		   &RTT_
		   ) //1.R
                                
		  + nut_*twoSymm(fvc::grad(U_))
        )
    )



//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//

{
    
    nut_ = Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    printCoeffs();	

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//

tmp<volTensorField> kEpsH::S_() const   
{
    return tmp<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "S_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           (Ltt_&symm(fvc::grad(U_)))&Rtt_
        )
    );
}

tmp<volTensorField> kEpsH::S() const
{
    return tmp<volTensorField>
    (
        new volTensorField
        (
            IOobject
            (
                "S",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           (Ltt_&symm(gradU))&Rtt_
        )
    );
}


tmp<volScalarField> kEpsH::SrPhi_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SrPhi_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::XY)
        )
    );
}

tmp<volScalarField> kEpsH::SrPhi() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SrPhi",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ), 
           S()->component(symmTensor::XY)
        )
    );
}


tmp<volScalarField> kEpsH::Srz_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Srz_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::XZ)
        )
    );
}

tmp<volScalarField> kEpsH::Srz() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Srz",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S()->component(symmTensor::XZ)
        )
    );
}


tmp<volScalarField> kEpsH::SPhiz_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SPhiz_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::YZ)
        )
    );
}

tmp<volScalarField> kEpsH::SPhiz() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SPhiz",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::YZ)
        )
    );
}


tmp<volScalarField> kEpsH::Szz_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Szz_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::ZZ)
        )
    );
}

tmp<volScalarField> kEpsH::Szz() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Szz",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S()->component(symmTensor::ZZ)
        )
    );
}

tmp<volScalarField> kEpsH::Srr_() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Srr_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::XX)
        )
    );
}

tmp<volScalarField> kEpsH::Srr() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Srr",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S()->component(symmTensor::XX)
        )
    );
}

tmp<volScalarField> kEpsH::SPhiPhi_() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SPhiPhi_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S_()->component(symmTensor::YY)
        )
    );
}

tmp<volScalarField> kEpsH::SPhiPhi() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "SPhiPhi",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           S()->component(symmTensor::YY)
        )
    );
}

tmp<volScalarField> kEpsH::Ux_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Ux_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           U_.component(vector::X)
        )
    );
}

tmp<volScalarField> kEpsH::Uy_() const    
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Uy_",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           U_.component(vector::Y)
        )
    );

}
//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +//

 tmp<volSymmTensorField> kEpsH::R() const
{  
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
	   IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + nonlinStress_,
            k_.boundaryField().types()
        )
    );

}

tmp<volSymmTensorField> kEpsH::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinStress_
        )
    );
}


tmp<fvVectorMatrix> kEpsH::divDevReff() const

{
    return
    (
              fvc::div(nonlinStress_)
            - fvm::laplacian(nuEff(), U_)
            - fvc::div(nuEff()*dev(T(fvc::grad(U_))))
    );
}



bool kEpsH::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsH::correct()
{

    if (mesh_.changing())
    {
        bound(k_, k0_);
        bound(epsilon_, epsilon0_);
    }

    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }
  
    gradU=fvc::grad(U_);  // DjN, 22/Aug/2018

    // generation term
    volScalarField G("RASModel::G", nut_*2*magSqr(symm(gradU)) - (nonlinStress_ && gradU));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();


    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    
    nut_ = Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_);
    nut_ = min(nut_, nuRatio()*nu());
    nut_.correctBoundaryConditions();

    //Recalculating Hamba's coefficients

    Omega_=mag(-Ux_()*SinPhi_+Uy_()*CosPhi_)/r_;
    TauH_=Ctau_*k_/(epsilon_ + epsilonSmall_);
    OT_=Omega_*TauH_;
    K1_=4*(OT_)/(1+4*sqr(OT_));
    K2_=2/(1+4*sqr(OT_));
    K3_=2/(1+sqr(OT_));
    K4_=2*(OT_)/(1+sqr(OT_));
    K5_=2*(OT_)/(1+4*sqr(OT_));
    K6_=4*sqr(OT_)/(1+4*sqr(OT_));
    K7_=2*(1+2*sqr(OT_))/(1+4*sqr(OT_));

    nonlinStress_ =
    symm(nut_*symm( //1.L
                   LTT_&
                   ( //2.L
                   ( //3.L
                   ((
                   -K1_*OneOne(tensor::one)
                   -K2_*OneTwo(tensor::one)
                   -K2_*TwoOne(tensor::one)
                   +K1_*TwoTwo(tensor::one)
                   ) 
                   &
                   (
                   SrPhi()*I
                   ))
                   +
                   ((
                   -K3_*OneThree(tensor::one)
                   +K4_*TwoThree(tensor::one)
                   -K3_*ThreeOne(tensor::one)
                   +K4_*ThreeTwo(tensor::one)
                   )
                   &
                   (
                   Srz()*I
                   ))
                   +
                   ((
                   -K4_*OneThree(tensor::one)
                   -K3_*TwoThree(tensor::one)
                   -K4_*ThreeOne(tensor::one)
                   -K3_*ThreeTwo(tensor::one)
                   )
                   &
                   (
                   SPhiz()*I
                   ))
                   +
                   ((
                   -2*ThreeThree(tensor::one)
                   )
                   &
                   (
                   Szz()*ThreeThree(tensor::one)
                   ))
		   +
		   ((
                   -K7_*OneOne(tensor::one)
                   +K5_*OneTwo(tensor::one)
                   +K5_*TwoOne(tensor::one)
                   -K6_*TwoTwo(tensor::one)
                   )
                   &
                   (
                   Srr()*I
		   ))
		   +
                   ((
                   -K6_*OneOne(tensor::one)
                   -K5_*OneTwo(tensor::one)
                   -K5_*TwoOne(tensor::one)
                   -K7_*TwoTwo(tensor::one)
                   )
                   &
                   (
                   SPhiPhi()*I
                   ))
                   ) //3.R
                   ) //2.R
                   &RTT_
                   ) //1.R
                             
                   + nut_*twoSymm(gradU)
    );


Info<<endl<<"You are using kEpsH model:"<<endl;
Info<<Ctau_<<endl<<endl;
Info<<Cmu_<<endl<<endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
