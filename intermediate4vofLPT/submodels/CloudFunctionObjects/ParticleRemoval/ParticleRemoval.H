/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::ParticleRemoval

Description
    Remove particles when they escape the alpha phase fraction; for multi-phase cases.

    Model is activated using:
    \verbatim
    particleRemoval1
    {
        type        particleRemoval;
        alphaName   alpha;      // name volume fraction field
        threshold   0.95;       // alpha value below which model is active
    }
    \endverbatim


SourceFiles
    ParticleRemoval.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleRemoval_H
#define ParticleRemoval_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleRemoval Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleRemoval
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Name of vol fraction field
        const word alphaName_;

        //- Pointer to the volume fraction field
        const volScalarField* alphaPtr_;

        //- Gradient of the volume fraction field
    // COMMENTED    autoPtr<volVectorField> gradAlphaPtr_;

        //- Threshold beyond which model is active
        scalar threshold_;


public:

    //- Runtime type information
    TypeName("particleRemoval");


    // Constructors

        //- Construct from dictionary
        ParticleRemoval
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleRemoval(const ParticleRemoval<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleRemoval<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleRemoval();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            virtual void preEvolveParticleRemoval
            (
                typename CloudType::parcelType& p,
                const label celli,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (
                typename CloudType::parcelType& p,
                const label celli,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

	    // ADDED LINES
            //- Post-patch hook
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const polyPatch& pp,
                const scalar trackFraction,
                const tetIndices& testIs,
                bool& keepParticle
            );

            //- Post-face hook
            virtual void postFace
            (
                const typename CloudType::parcelType& p,
                const label faceI,
                bool& keepParticle
            );
            // END ADDED LINES	    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleRemoval.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
