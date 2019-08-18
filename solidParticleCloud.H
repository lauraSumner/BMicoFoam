/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::solidParticleCloud

Description
    A Cloud of solid particles

SourceFiles
    solidParticleCloudI.H
    solidParticleCloud.C
    solidParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCloud_H
#define solidParticleCloud_H

#include "Cloud.H"
#include "solidParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class solidParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCloud
:
    public Cloud<solidParticle>
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        scalar rhop_;
        scalar e_;
        scalar mu_;
		vector posP1_;
		scalar dP1_;
		vector UP1_;
		vector posP2_;
		scalar dP2_;
		vector UP2_;
		vector posP3_;
		scalar dP3_;
		vector UP3_;
		vector posP4_;
		scalar dP4_;
		vector UP4_;
		vector posP5_;
		scalar dP5_;
		vector UP5_;
		vector posP6_;
		scalar dP6_;
		vector UP6_;
		vector posP7_;
		scalar dP7_;
		vector UP7_;
		vector posP8_;
		scalar dP8_;
		vector UP8_;
		vector posP9_;
		scalar dP9_;
		vector UP9_;
		vector posP10_;
		scalar dP10_;
		vector UP10_;
		vector posP11_;
		scalar dP11_;
		vector UP11_;
		vector posP12_;
		scalar dP12_;
		vector UP12_;
		vector posP13_;
		scalar dP13_;
		vector UP13_;
		vector posP14_;
		scalar dP14_;
		vector UP14_;
		vector posP15_;
		scalar dP15_;
		vector UP15_;
		vector posP16_;
		scalar dP16_;
		vector UP16_;
		vector posP17_;
		scalar dP17_;
		vector UP17_;
		vector posP18_;
		scalar dP18_;
		vector UP18_;
		vector posP19_;
		scalar dP19_;
		vector UP19_;
		vector posP20_;
		scalar dP20_;
		vector UP20_;
		vector posP21_;
		scalar dP21_;
		vector UP21_;
		vector posP22_;
		scalar dP22_;
		vector UP22_;
		vector posP23_;
		scalar dP23_;
		vector UP23_;
		vector posP24_;
		scalar dP24_;
		vector UP24_;
		vector posP25_;
		scalar dP25_;
		vector UP25_;
		vector posP26_;
		scalar dP26_;
		vector UP26_;
		vector posP27_;
		scalar dP27_;
		vector UP27_;
		vector posP28_;
		scalar dP28_;
		vector UP28_;
		vector posP29_;
		scalar dP29_;
		vector UP29_;
		vector posP30_;
		scalar dP30_;
		vector UP30_;
		vector posP31_;
		scalar dP31_;
		vector UP31_;
		vector posP32_;
		scalar dP32_;
		vector UP32_;
		vector posP33_;
		scalar dP33_;
		vector UP33_;
		vector posP34_;
		scalar dP34_;
		vector UP34_;
		vector posP35_;
		scalar dP35_;
		vector UP35_;
		vector posP36_;
		scalar dP36_;
		vector UP36_;
		vector posP37_;
		scalar dP37_;
		vector UP37_;
		vector posP38_;
		scalar dP38_;
		vector UP38_;
		vector posP39_;
		scalar dP39_;
		vector UP39_;
		vector posP40_;
		scalar dP40_;
		vector UP40_;
		vector posP41_;
		scalar dP41_;
		vector UP41_;
		vector posP42_;
		scalar dP42_;
		vector UP42_;
		vector posP43_;
		scalar dP43_;
		vector UP43_;
		vector posP44_;
		scalar dP44_;
		vector UP44_;
		vector posP45_;
		scalar dP45_;
		vector UP45_;
		vector posP46_;
		scalar dP46_;
		vector UP46_;
		vector posP47_;
		scalar dP47_;
		vector UP47_;
		vector posP48_;
		scalar dP48_;
		vector UP48_;
		vector posP49_;
		scalar dP49_;
		vector UP49_;
		vector posP50_;
		scalar dP50_;
		vector UP50_;
		vector posP51_;
		scalar dP51_;
		vector UP51_;
		vector posP52_;
		scalar dP52_;
		vector UP52_;
		vector posP53_;
		scalar dP53_;
		vector UP53_;
		vector posP54_;
		scalar dP54_;
		vector UP54_;
		vector posP55_;
		scalar dP55_;
		vector UP55_;
		vector posP56_;
		scalar dP56_;
		vector UP56_;
		vector posP57_;
		scalar dP57_;
		vector UP57_;
		vector posP58_;
		scalar dP58_;
		vector UP58_;
		vector posP59_;
		scalar dP59_;
		vector UP59_;
		vector posP60_;
		scalar dP60_;
		vector UP60_;
		vector posP61_;
		scalar dP61_;
		vector UP61_;
		vector posP62_;
		scalar dP62_;
		vector UP62_;
		vector posP63_;
		scalar dP63_;
		vector UP63_;
		vector posP64_;
		scalar dP64_;
		vector UP64_;
		vector posP65_;
		scalar dP65_;
		vector UP65_;
		vector posP66_;
		scalar dP66_;
		vector UP66_;
		vector posP67_;
		scalar dP67_;
		vector UP67_;
		vector posP68_;
		scalar dP68_;
		vector UP68_;
		vector posP69_;
		scalar dP69_;
		vector UP69_;
		vector posP70_;
		scalar dP70_;
		vector UP70_;
		vector posP71_;
		scalar dP71_;
		vector UP71_;
		vector posP72_;
		scalar dP72_;
		vector UP72_;
		vector posP73_;
		scalar dP73_;
		vector UP73_;
		vector posP74_;
		scalar dP74_;
		vector UP74_;
		vector posP75_;
		scalar dP75_;
		vector UP75_;
		vector posP76_;
		scalar dP76_;
		vector UP76_;
		vector posP77_;
		scalar dP77_;
		vector UP77_;
		vector posP78_;
		scalar dP78_;
		vector UP78_;
		vector posP79_;
		scalar dP79_;
		vector UP79_;
		vector posP80_;
		scalar dP80_;
		vector UP80_;
		vector posP81_;
		scalar dP81_;
		vector UP81_;
		vector posP82_;
		scalar dP82_;
		vector UP82_;
		vector posP83_;
		scalar dP83_;
		vector UP83_;
		vector posP84_;
		scalar dP84_;
		vector UP84_;
		vector posP85_;
		scalar dP85_;
		vector UP85_;
		vector posP86_;
		scalar dP86_;
		vector UP86_;
		vector posP87_;
		scalar dP87_;
		vector UP87_;
		vector posP88_;
		scalar dP88_;
		vector UP88_;
		vector posP89_;
		scalar dP89_;
		vector UP89_;
		vector posP90_;
		scalar dP90_;
		vector UP90_;
		vector posP91_;
		scalar dP91_;
		vector UP91_;
		vector posP92_;
		scalar dP92_;
		vector UP92_;
		vector posP93_;
		scalar dP93_;
		vector UP93_;
		vector posP94_;
		scalar dP94_;
		vector UP94_;
		vector posP95_;
		scalar dP95_;
		vector UP95_;
		vector posP96_;
		scalar dP96_;
		vector UP96_;
		vector posP97_;
		scalar dP97_;
		vector UP97_;
		vector posP98_;
		scalar dP98_;
		vector UP98_;
		vector posP99_;
		scalar dP99_;
		vector UP99_;
		vector posP100_;
		scalar dP100_;
		vector UP100_;
		
		

		scalar tInjStart_;
		scalar tInjEnd_;
		


    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidParticleCloud(const solidParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const solidParticleCloud&);


public:

    // Constructors

        //- Construct given mesh
        solidParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );


    // Member Functions

        // Access

            virtual bool hasWallImpactDistance() const;

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mu() const;


        // Edit

            //- Move the particles under the influence of the given
            //  gravitational acceleration
            void move(const dimensionedVector& g);
//- Inject particles according to the dictionnary particleProperties
			void inject(solidParticle::trackingData &td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
