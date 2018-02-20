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

\*---------------------------------------------------------------------------*/

#include "ParticleRemoval.H"
#include "fvcGrad.H"

#define PARTICLEREMOVAL_DEBUG true

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRemoval<CloudType>::ParticleRemoval
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    alphaName_
    (
        this->coeffDict().template lookupOrDefault<word>("alpha", "alpha")
    ),
    alphaPtr_(NULL),
    // COMMENTED gradAlphaPtr_(NULL),
    threshold_(readScalar(this->coeffDict().lookup("threshold")))
{}


template<class CloudType>
Foam::ParticleRemoval<CloudType>::ParticleRemoval
(
    const ParticleRemoval<CloudType>& pt
)
:
    CloudFunctionObject<CloudType>(pt),
    alphaName_(pt.alphaName_),
    alphaPtr_(pt.alphaPtr_),
    // COMMENTED gradAlphaPtr_(NULL),
    threshold_(pt.threshold_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRemoval<CloudType>::~ParticleRemoval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleRemoval<CloudType>::preEvolve()
{
    if (alphaPtr_ == NULL)
    {
        const fvMesh& mesh = this->owner().mesh();
        const volScalarField& alpha =
            mesh.lookupObject<volScalarField>(alphaName_);

        alphaPtr_ = &alpha;
    }
}

template<class CloudType>
void Foam::ParticleRemoval<CloudType>::preEvolveParticleRemoval
(
    parcelType& p,
    const label celli,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    bool& active = p.active();

    if (active == false && keepParticle == true)
    {
        keepParticle = false;
#if PARTICLEREMOVAL_DEBUG
        Pout<< "ParticleRemoval-keepParticle is = " << keepParticle << nl << endl;
        Pout<< "ParticleRemoval-active is = " << active << nl << endl;
#endif
    }
}


template<class CloudType>
void Foam::ParticleRemoval<CloudType>::postEvolve()
{
    // do nothing
}


template<class CloudType>
void Foam::ParticleRemoval<CloudType>::postMove
(
    parcelType& p,
    const label celli,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    bool& active = p.active();

    if (alphaPtr_->primitiveField()[celli] < threshold_ && active == true)
    {
        keepParticle = true; //false;
        active = false;
#if PARTICLEREMOVAL_DEBUG
        Pout<< "ParticleRemoval-keepParticle2 is = " << keepParticle << nl << endl;
        Pout<< "ParticleRemoval-active2 is = " << active << nl << endl;
#endif
    }
}

// ADDED LINES
template<class CloudType>
void Foam::ParticleRemoval<CloudType>::postPatch
(
    const typename CloudType::parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& testIs,
    bool& keepParticle
)
{
    // Do nothing
}


template<class CloudType>
void Foam::ParticleRemoval<CloudType>::postFace
(
    const typename CloudType::parcelType& p,
    const label faceI,
    bool& keepParticle
)
{
    // Do nothing
}
// END ADDED LINES

// ************************************************************************* //
