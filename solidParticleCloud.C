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

\*---------------------------------------------------------------------------*/

#include "solidParticleCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticleCloud::solidParticleCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    bool readFields
)
:
    Cloud<solidParticle>(mesh, cloudName, false),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            "particleProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    rhop_(dimensionedScalar(particleProperties_.lookup("rhop")).value()),
    e_(dimensionedScalar(particleProperties_.lookup("e")).value()),
    mu_(dimensionedScalar(particleProperties_.lookup("mu")).value()),
	posP1_(dimensionedVector(particleProperties_.lookup("posP1")).value()),
	dP1_(dimensionedScalar(particleProperties_.lookup("dP1")).value()),
	UP1_(dimensionedVector(particleProperties_.lookup("UP1")).value()),
	posP2_(dimensionedVector(particleProperties_.lookup("posP2")).value()),
	dP2_(dimensionedScalar(particleProperties_.lookup("dP2")).value()),
	UP2_(dimensionedVector(particleProperties_.lookup("UP2")).value()),
	posP3_(dimensionedVector(particleProperties_.lookup("posP3")).value()),
	dP3_(dimensionedScalar(particleProperties_.lookup("dP3")).value()),
	UP3_(dimensionedVector(particleProperties_.lookup("UP3")).value()),
	posP4_(dimensionedVector(particleProperties_.lookup("posP4")).value()),
	dP4_(dimensionedScalar(particleProperties_.lookup("dP4")).value()),
	UP4_(dimensionedVector(particleProperties_.lookup("UP4")).value()),
	posP5_(dimensionedVector(particleProperties_.lookup("posP5")).value()),
	dP5_(dimensionedScalar(particleProperties_.lookup("dP5")).value()),
	UP5_(dimensionedVector(particleProperties_.lookup("UP5")).value()),
	posP6_(dimensionedVector(particleProperties_.lookup("posP6")).value()),
	dP6_(dimensionedScalar(particleProperties_.lookup("dP6")).value()),
	UP6_(dimensionedVector(particleProperties_.lookup("UP6")).value()),
	posP7_(dimensionedVector(particleProperties_.lookup("posP7")).value()),
	dP7_(dimensionedScalar(particleProperties_.lookup("dP7")).value()),
	UP7_(dimensionedVector(particleProperties_.lookup("UP7")).value()),
	posP8_(dimensionedVector(particleProperties_.lookup("posP8")).value()),
	dP8_(dimensionedScalar(particleProperties_.lookup("dP8")).value()),
	UP8_(dimensionedVector(particleProperties_.lookup("UP8")).value()),
	posP9_(dimensionedVector(particleProperties_.lookup("posP9")).value()),
	dP9_(dimensionedScalar(particleProperties_.lookup("dP9")).value()),
	UP9_(dimensionedVector(particleProperties_.lookup("UP9")).value()),
	posP10_(dimensionedVector(particleProperties_.lookup("posP10")).value()),
	dP10_(dimensionedScalar(particleProperties_.lookup("dP10")).value()),
	UP10_(dimensionedVector(particleProperties_.lookup("UP10")).value()),
	posP11_(dimensionedVector(particleProperties_.lookup("posP11")).value()),
	dP11_(dimensionedScalar(particleProperties_.lookup("dP11")).value()),
	UP11_(dimensionedVector(particleProperties_.lookup("UP11")).value()),
	posP12_(dimensionedVector(particleProperties_.lookup("posP12")).value()),
	dP12_(dimensionedScalar(particleProperties_.lookup("dP12")).value()),
	UP12_(dimensionedVector(particleProperties_.lookup("UP12")).value()),
	posP13_(dimensionedVector(particleProperties_.lookup("posP13")).value()),
	dP13_(dimensionedScalar(particleProperties_.lookup("dP13")).value()),
	UP13_(dimensionedVector(particleProperties_.lookup("UP13")).value()),
	posP14_(dimensionedVector(particleProperties_.lookup("posP14")).value()),
	dP14_(dimensionedScalar(particleProperties_.lookup("dP14")).value()),
	UP14_(dimensionedVector(particleProperties_.lookup("UP14")).value()),
	posP15_(dimensionedVector(particleProperties_.lookup("posP15")).value()),
	dP15_(dimensionedScalar(particleProperties_.lookup("dP15")).value()),
	UP15_(dimensionedVector(particleProperties_.lookup("UP15")).value()),
	posP16_(dimensionedVector(particleProperties_.lookup("posP16")).value()),
	dP16_(dimensionedScalar(particleProperties_.lookup("dP16")).value()),
	UP16_(dimensionedVector(particleProperties_.lookup("UP16")).value()),
	posP17_(dimensionedVector(particleProperties_.lookup("posP17")).value()),
	dP17_(dimensionedScalar(particleProperties_.lookup("dP17")).value()),
	UP17_(dimensionedVector(particleProperties_.lookup("UP17")).value()),
	posP18_(dimensionedVector(particleProperties_.lookup("posP18")).value()),
	dP18_(dimensionedScalar(particleProperties_.lookup("dP18")).value()),
	UP18_(dimensionedVector(particleProperties_.lookup("UP18")).value()),
	posP19_(dimensionedVector(particleProperties_.lookup("posP19")).value()),
	dP19_(dimensionedScalar(particleProperties_.lookup("dP19")).value()),
	UP19_(dimensionedVector(particleProperties_.lookup("UP19")).value()),
	posP20_(dimensionedVector(particleProperties_.lookup("posP20")).value()),
	dP20_(dimensionedScalar(particleProperties_.lookup("dP20")).value()),
	UP20_(dimensionedVector(particleProperties_.lookup("UP20")).value()),
	posP21_(dimensionedVector(particleProperties_.lookup("posP21")).value()),
	dP21_(dimensionedScalar(particleProperties_.lookup("dP21")).value()),
	UP21_(dimensionedVector(particleProperties_.lookup("UP21")).value()),
	posP22_(dimensionedVector(particleProperties_.lookup("posP22")).value()),
	dP22_(dimensionedScalar(particleProperties_.lookup("dP22")).value()),
	UP22_(dimensionedVector(particleProperties_.lookup("UP22")).value()),
	posP23_(dimensionedVector(particleProperties_.lookup("posP23")).value()),
	dP23_(dimensionedScalar(particleProperties_.lookup("dP23")).value()),
	UP23_(dimensionedVector(particleProperties_.lookup("UP23")).value()),
	posP24_(dimensionedVector(particleProperties_.lookup("posP24")).value()),
	dP24_(dimensionedScalar(particleProperties_.lookup("dP24")).value()),
	UP24_(dimensionedVector(particleProperties_.lookup("UP24")).value()),
	posP25_(dimensionedVector(particleProperties_.lookup("posP25")).value()),
	dP25_(dimensionedScalar(particleProperties_.lookup("dP25")).value()),
	UP25_(dimensionedVector(particleProperties_.lookup("UP25")).value()),
	posP26_(dimensionedVector(particleProperties_.lookup("posP26")).value()),
	dP26_(dimensionedScalar(particleProperties_.lookup("dP26")).value()),
	UP26_(dimensionedVector(particleProperties_.lookup("UP26")).value()),
	posP27_(dimensionedVector(particleProperties_.lookup("posP27")).value()),
	dP27_(dimensionedScalar(particleProperties_.lookup("dP27")).value()),
	UP27_(dimensionedVector(particleProperties_.lookup("UP27")).value()),
	posP28_(dimensionedVector(particleProperties_.lookup("posP28")).value()),
	dP28_(dimensionedScalar(particleProperties_.lookup("dP28")).value()),
	UP28_(dimensionedVector(particleProperties_.lookup("UP28")).value()),
	posP29_(dimensionedVector(particleProperties_.lookup("posP29")).value()),
	dP29_(dimensionedScalar(particleProperties_.lookup("dP29")).value()),
	UP29_(dimensionedVector(particleProperties_.lookup("UP29")).value()),
	posP30_(dimensionedVector(particleProperties_.lookup("posP30")).value()),
	dP30_(dimensionedScalar(particleProperties_.lookup("dP30")).value()),
	UP30_(dimensionedVector(particleProperties_.lookup("UP30")).value()),
	posP31_(dimensionedVector(particleProperties_.lookup("posP31")).value()),
	dP31_(dimensionedScalar(particleProperties_.lookup("dP31")).value()),
	UP31_(dimensionedVector(particleProperties_.lookup("UP31")).value()),
	posP32_(dimensionedVector(particleProperties_.lookup("posP32")).value()),
	dP32_(dimensionedScalar(particleProperties_.lookup("dP32")).value()),
	UP32_(dimensionedVector(particleProperties_.lookup("UP32")).value()),
	posP33_(dimensionedVector(particleProperties_.lookup("posP33")).value()),
	dP33_(dimensionedScalar(particleProperties_.lookup("dP33")).value()),
	UP33_(dimensionedVector(particleProperties_.lookup("UP33")).value()),
	posP34_(dimensionedVector(particleProperties_.lookup("posP34")).value()),
	dP34_(dimensionedScalar(particleProperties_.lookup("dP34")).value()),
	UP34_(dimensionedVector(particleProperties_.lookup("UP34")).value()),
	posP35_(dimensionedVector(particleProperties_.lookup("posP35")).value()),
	dP35_(dimensionedScalar(particleProperties_.lookup("dP35")).value()),
	UP35_(dimensionedVector(particleProperties_.lookup("UP35")).value()),
	posP36_(dimensionedVector(particleProperties_.lookup("posP36")).value()),
	dP36_(dimensionedScalar(particleProperties_.lookup("dP36")).value()),
	UP36_(dimensionedVector(particleProperties_.lookup("UP36")).value()),
	posP37_(dimensionedVector(particleProperties_.lookup("posP37")).value()),
	dP37_(dimensionedScalar(particleProperties_.lookup("dP37")).value()),
	UP37_(dimensionedVector(particleProperties_.lookup("UP37")).value()),
	posP38_(dimensionedVector(particleProperties_.lookup("posP38")).value()),
	dP38_(dimensionedScalar(particleProperties_.lookup("dP38")).value()),
	UP38_(dimensionedVector(particleProperties_.lookup("UP38")).value()),
	posP39_(dimensionedVector(particleProperties_.lookup("posP39")).value()),
	dP39_(dimensionedScalar(particleProperties_.lookup("dP39")).value()),
	UP39_(dimensionedVector(particleProperties_.lookup("UP39")).value()),
	posP40_(dimensionedVector(particleProperties_.lookup("posP40")).value()),
	dP40_(dimensionedScalar(particleProperties_.lookup("dP40")).value()),
	UP40_(dimensionedVector(particleProperties_.lookup("UP40")).value()),
	posP41_(dimensionedVector(particleProperties_.lookup("posP41")).value()),
	dP41_(dimensionedScalar(particleProperties_.lookup("dP41")).value()),
	UP41_(dimensionedVector(particleProperties_.lookup("UP41")).value()),
	posP42_(dimensionedVector(particleProperties_.lookup("posP42")).value()),
	dP42_(dimensionedScalar(particleProperties_.lookup("dP42")).value()),
	UP42_(dimensionedVector(particleProperties_.lookup("UP42")).value()),
	posP43_(dimensionedVector(particleProperties_.lookup("posP43")).value()),
	dP43_(dimensionedScalar(particleProperties_.lookup("dP43")).value()),
	UP43_(dimensionedVector(particleProperties_.lookup("UP43")).value()),
	posP44_(dimensionedVector(particleProperties_.lookup("posP44")).value()),
	dP44_(dimensionedScalar(particleProperties_.lookup("dP44")).value()),
	UP44_(dimensionedVector(particleProperties_.lookup("UP44")).value()),
	posP45_(dimensionedVector(particleProperties_.lookup("posP45")).value()),
	dP45_(dimensionedScalar(particleProperties_.lookup("dP45")).value()),
	UP45_(dimensionedVector(particleProperties_.lookup("UP45")).value()),
	posP46_(dimensionedVector(particleProperties_.lookup("posP46")).value()),
	dP46_(dimensionedScalar(particleProperties_.lookup("dP46")).value()),
	UP46_(dimensionedVector(particleProperties_.lookup("UP46")).value()),
	posP47_(dimensionedVector(particleProperties_.lookup("posP47")).value()),
	dP47_(dimensionedScalar(particleProperties_.lookup("dP47")).value()),
	UP47_(dimensionedVector(particleProperties_.lookup("UP47")).value()),
	posP48_(dimensionedVector(particleProperties_.lookup("posP48")).value()),
	dP48_(dimensionedScalar(particleProperties_.lookup("dP48")).value()),
	UP48_(dimensionedVector(particleProperties_.lookup("UP48")).value()),
	posP49_(dimensionedVector(particleProperties_.lookup("posP49")).value()),
	dP49_(dimensionedScalar(particleProperties_.lookup("dP49")).value()),
	UP49_(dimensionedVector(particleProperties_.lookup("UP49")).value()),
	posP50_(dimensionedVector(particleProperties_.lookup("posP50")).value()),
	dP50_(dimensionedScalar(particleProperties_.lookup("dP50")).value()),
	UP50_(dimensionedVector(particleProperties_.lookup("UP50")).value()),
	posP51_(dimensionedVector(particleProperties_.lookup("posP51")).value()),
	dP51_(dimensionedScalar(particleProperties_.lookup("dP51")).value()),
	UP51_(dimensionedVector(particleProperties_.lookup("UP51")).value()),
	posP52_(dimensionedVector(particleProperties_.lookup("posP52")).value()),
	dP52_(dimensionedScalar(particleProperties_.lookup("dP52")).value()),
	UP52_(dimensionedVector(particleProperties_.lookup("UP52")).value()),
	posP53_(dimensionedVector(particleProperties_.lookup("posP53")).value()),
	dP53_(dimensionedScalar(particleProperties_.lookup("dP53")).value()),
	UP53_(dimensionedVector(particleProperties_.lookup("UP53")).value()),
	posP54_(dimensionedVector(particleProperties_.lookup("posP54")).value()),
	dP54_(dimensionedScalar(particleProperties_.lookup("dP54")).value()),
	UP54_(dimensionedVector(particleProperties_.lookup("UP54")).value()),
	posP55_(dimensionedVector(particleProperties_.lookup("posP55")).value()),
	dP55_(dimensionedScalar(particleProperties_.lookup("dP55")).value()),
	UP55_(dimensionedVector(particleProperties_.lookup("UP55")).value()),
	posP56_(dimensionedVector(particleProperties_.lookup("posP56")).value()),
	dP56_(dimensionedScalar(particleProperties_.lookup("dP56")).value()),
	UP56_(dimensionedVector(particleProperties_.lookup("UP56")).value()),
	posP57_(dimensionedVector(particleProperties_.lookup("posP57")).value()),
	dP57_(dimensionedScalar(particleProperties_.lookup("dP57")).value()),
	UP57_(dimensionedVector(particleProperties_.lookup("UP57")).value()),
	posP58_(dimensionedVector(particleProperties_.lookup("posP58")).value()),
	dP58_(dimensionedScalar(particleProperties_.lookup("dP58")).value()),
	UP58_(dimensionedVector(particleProperties_.lookup("UP58")).value()),
	posP59_(dimensionedVector(particleProperties_.lookup("posP59")).value()),
	dP59_(dimensionedScalar(particleProperties_.lookup("dP59")).value()),
	UP59_(dimensionedVector(particleProperties_.lookup("UP59")).value()),
	posP60_(dimensionedVector(particleProperties_.lookup("posP60")).value()),
	dP60_(dimensionedScalar(particleProperties_.lookup("dP60")).value()),
	UP60_(dimensionedVector(particleProperties_.lookup("UP60")).value()),
	posP61_(dimensionedVector(particleProperties_.lookup("posP61")).value()),
	dP61_(dimensionedScalar(particleProperties_.lookup("dP61")).value()),
	UP61_(dimensionedVector(particleProperties_.lookup("UP61")).value()),
	posP62_(dimensionedVector(particleProperties_.lookup("posP62")).value()),
	dP62_(dimensionedScalar(particleProperties_.lookup("dP62")).value()),
	UP62_(dimensionedVector(particleProperties_.lookup("UP62")).value()),
	posP63_(dimensionedVector(particleProperties_.lookup("posP63")).value()),
	dP63_(dimensionedScalar(particleProperties_.lookup("dP63")).value()),
	UP63_(dimensionedVector(particleProperties_.lookup("UP63")).value()),
	posP64_(dimensionedVector(particleProperties_.lookup("posP64")).value()),
	dP64_(dimensionedScalar(particleProperties_.lookup("dP64")).value()),
	UP64_(dimensionedVector(particleProperties_.lookup("UP64")).value()),
	posP65_(dimensionedVector(particleProperties_.lookup("posP65")).value()),
	dP65_(dimensionedScalar(particleProperties_.lookup("dP65")).value()),
	UP65_(dimensionedVector(particleProperties_.lookup("UP65")).value()),
	posP66_(dimensionedVector(particleProperties_.lookup("posP66")).value()),
	dP66_(dimensionedScalar(particleProperties_.lookup("dP66")).value()),
	UP66_(dimensionedVector(particleProperties_.lookup("UP66")).value()),
	posP67_(dimensionedVector(particleProperties_.lookup("posP67")).value()),
	dP67_(dimensionedScalar(particleProperties_.lookup("dP67")).value()),
	UP67_(dimensionedVector(particleProperties_.lookup("UP67")).value()),
	posP68_(dimensionedVector(particleProperties_.lookup("posP68")).value()),
	dP68_(dimensionedScalar(particleProperties_.lookup("dP68")).value()),
	UP68_(dimensionedVector(particleProperties_.lookup("UP68")).value()),
	posP69_(dimensionedVector(particleProperties_.lookup("posP69")).value()),
	dP69_(dimensionedScalar(particleProperties_.lookup("dP69")).value()),
	UP69_(dimensionedVector(particleProperties_.lookup("UP69")).value()),
	posP70_(dimensionedVector(particleProperties_.lookup("posP70")).value()),
	dP70_(dimensionedScalar(particleProperties_.lookup("dP70")).value()),
	UP70_(dimensionedVector(particleProperties_.lookup("UP70")).value()),
	posP71_(dimensionedVector(particleProperties_.lookup("posP71")).value()),
	dP71_(dimensionedScalar(particleProperties_.lookup("dP71")).value()),
	UP71_(dimensionedVector(particleProperties_.lookup("UP71")).value()),
	posP72_(dimensionedVector(particleProperties_.lookup("posP72")).value()),
	dP72_(dimensionedScalar(particleProperties_.lookup("dP72")).value()),
	UP72_(dimensionedVector(particleProperties_.lookup("UP72")).value()),
	posP73_(dimensionedVector(particleProperties_.lookup("posP73")).value()),
	dP73_(dimensionedScalar(particleProperties_.lookup("dP73")).value()),
	UP73_(dimensionedVector(particleProperties_.lookup("UP73")).value()),
	posP74_(dimensionedVector(particleProperties_.lookup("posP74")).value()),
	dP74_(dimensionedScalar(particleProperties_.lookup("dP74")).value()),
	UP74_(dimensionedVector(particleProperties_.lookup("UP74")).value()),
	posP75_(dimensionedVector(particleProperties_.lookup("posP75")).value()),
	dP75_(dimensionedScalar(particleProperties_.lookup("dP75")).value()),
	UP75_(dimensionedVector(particleProperties_.lookup("UP75")).value()),
	posP76_(dimensionedVector(particleProperties_.lookup("posP76")).value()),
	dP76_(dimensionedScalar(particleProperties_.lookup("dP76")).value()),
	UP76_(dimensionedVector(particleProperties_.lookup("UP76")).value()),
	posP77_(dimensionedVector(particleProperties_.lookup("posP77")).value()),
	dP77_(dimensionedScalar(particleProperties_.lookup("dP77")).value()),
	UP77_(dimensionedVector(particleProperties_.lookup("UP77")).value()),
	posP78_(dimensionedVector(particleProperties_.lookup("posP78")).value()),
	dP78_(dimensionedScalar(particleProperties_.lookup("dP78")).value()),
	UP78_(dimensionedVector(particleProperties_.lookup("UP78")).value()),
	posP79_(dimensionedVector(particleProperties_.lookup("posP79")).value()),
	dP79_(dimensionedScalar(particleProperties_.lookup("dP79")).value()),
	UP79_(dimensionedVector(particleProperties_.lookup("UP79")).value()),
	posP80_(dimensionedVector(particleProperties_.lookup("posP80")).value()),
	dP80_(dimensionedScalar(particleProperties_.lookup("dP80")).value()),
	UP80_(dimensionedVector(particleProperties_.lookup("UP80")).value()),
	posP81_(dimensionedVector(particleProperties_.lookup("posP81")).value()),
	dP81_(dimensionedScalar(particleProperties_.lookup("dP81")).value()),
	UP81_(dimensionedVector(particleProperties_.lookup("UP81")).value()),
	posP82_(dimensionedVector(particleProperties_.lookup("posP82")).value()),
	dP82_(dimensionedScalar(particleProperties_.lookup("dP82")).value()),
	UP82_(dimensionedVector(particleProperties_.lookup("UP82")).value()),
	posP83_(dimensionedVector(particleProperties_.lookup("posP83")).value()),
	dP83_(dimensionedScalar(particleProperties_.lookup("dP83")).value()),
	UP83_(dimensionedVector(particleProperties_.lookup("UP83")).value()),
	posP84_(dimensionedVector(particleProperties_.lookup("posP84")).value()),
	dP84_(dimensionedScalar(particleProperties_.lookup("dP84")).value()),
	UP84_(dimensionedVector(particleProperties_.lookup("UP84")).value()),
	posP85_(dimensionedVector(particleProperties_.lookup("posP85")).value()),
	dP85_(dimensionedScalar(particleProperties_.lookup("dP85")).value()),
	UP85_(dimensionedVector(particleProperties_.lookup("UP85")).value()),
	posP86_(dimensionedVector(particleProperties_.lookup("posP86")).value()),
	dP86_(dimensionedScalar(particleProperties_.lookup("dP86")).value()),
	UP86_(dimensionedVector(particleProperties_.lookup("UP86")).value()),
	posP87_(dimensionedVector(particleProperties_.lookup("posP87")).value()),
	dP87_(dimensionedScalar(particleProperties_.lookup("dP87")).value()),
	UP87_(dimensionedVector(particleProperties_.lookup("UP87")).value()),
	posP88_(dimensionedVector(particleProperties_.lookup("posP88")).value()),
	dP88_(dimensionedScalar(particleProperties_.lookup("dP88")).value()),
	UP88_(dimensionedVector(particleProperties_.lookup("UP88")).value()),
	posP89_(dimensionedVector(particleProperties_.lookup("posP89")).value()),
	dP89_(dimensionedScalar(particleProperties_.lookup("dP89")).value()),
	UP89_(dimensionedVector(particleProperties_.lookup("UP89")).value()),
	posP90_(dimensionedVector(particleProperties_.lookup("posP90")).value()),
	dP90_(dimensionedScalar(particleProperties_.lookup("dP90")).value()),
	UP90_(dimensionedVector(particleProperties_.lookup("UP90")).value()),
	posP91_(dimensionedVector(particleProperties_.lookup("posP91")).value()),
	dP91_(dimensionedScalar(particleProperties_.lookup("dP91")).value()),
	UP91_(dimensionedVector(particleProperties_.lookup("UP91")).value()),
	posP92_(dimensionedVector(particleProperties_.lookup("posP92")).value()),
	dP92_(dimensionedScalar(particleProperties_.lookup("dP92")).value()),
	UP92_(dimensionedVector(particleProperties_.lookup("UP92")).value()),
	posP93_(dimensionedVector(particleProperties_.lookup("posP93")).value()),
	dP93_(dimensionedScalar(particleProperties_.lookup("dP93")).value()),
	UP93_(dimensionedVector(particleProperties_.lookup("UP93")).value()),
	posP94_(dimensionedVector(particleProperties_.lookup("posP94")).value()),
	dP94_(dimensionedScalar(particleProperties_.lookup("dP94")).value()),
	UP94_(dimensionedVector(particleProperties_.lookup("UP94")).value()),
	posP95_(dimensionedVector(particleProperties_.lookup("posP95")).value()),
	dP95_(dimensionedScalar(particleProperties_.lookup("dP95")).value()),
	UP95_(dimensionedVector(particleProperties_.lookup("UP95")).value()),
	posP96_(dimensionedVector(particleProperties_.lookup("posP96")).value()),
	dP96_(dimensionedScalar(particleProperties_.lookup("dP96")).value()),
	UP96_(dimensionedVector(particleProperties_.lookup("UP96")).value()),
	posP97_(dimensionedVector(particleProperties_.lookup("posP97")).value()),
	dP97_(dimensionedScalar(particleProperties_.lookup("dP97")).value()),
	UP97_(dimensionedVector(particleProperties_.lookup("UP97")).value()),
	posP98_(dimensionedVector(particleProperties_.lookup("posP98")).value()),
	dP98_(dimensionedScalar(particleProperties_.lookup("dP98")).value()),
	UP98_(dimensionedVector(particleProperties_.lookup("UP98")).value()),
	posP99_(dimensionedVector(particleProperties_.lookup("posP99")).value()),
	dP99_(dimensionedScalar(particleProperties_.lookup("dP99")).value()),
	UP99_(dimensionedVector(particleProperties_.lookup("UP99")).value()),
	posP100_(dimensionedVector(particleProperties_.lookup("posP100")).value()),
	dP100_(dimensionedScalar(particleProperties_.lookup("dP100")).value()),
	UP100_(dimensionedVector(particleProperties_.lookup("UP100")).value()),

	tInjStart_(dimensionedScalar(particleProperties_.lookup("tInjStart")).value()),
	tInjEnd_(dimensionedScalar(particleProperties_.lookup("tInjEnd")).value())

{	
    if (readFields)
    {
        solidParticle::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void Foam::solidParticleCloud::move(const dimensionedVector& g)
{
    const volScalarField& rho = mesh_.lookupObject<const volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
    const volScalarField& nu = mesh_.lookupObject<const volScalarField>("nu");

    interpolationCellPoint<scalar> rhoInterp(rho);
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> nuInterp(nu);
	
    solidParticle::trackingData
        td(*this, rhoInterp, UInterp, nuInterp, g.value());

    Cloud<solidParticle>::move(td, mesh_.time().deltaTValue());				
		if(mesh_.time().value()> td.cloud().tInjStart_ &&mesh_.time().value()< td.cloud().tInjEnd_ )
		{
		this->inject(td);
		}
}

void Foam::solidParticleCloud::inject(solidParticle::trackingData &td)
{
  label cellI=1;
  label tetFaceI=1;
  label tetPtI=1;

double  injectionX1coordinate = td.cloud().posP1_[0];	
double  injectionX2coordinate = td.cloud().posP2_[0];
double  injectionY1coordinate = td.cloud().posP1_[1];	
double  injectionY2coordinate = td.cloud().posP2_[1];
double  injectionX3coordinate = td.cloud().posP3_[0];
double  injectionX4coordinate = td.cloud().posP4_[0];
double  injectionX5coordinate = td.cloud().posP5_[0];
double  injectionX6coordinate = td.cloud().posP6_[0];
double  injectionX7coordinate = td.cloud().posP7_[0];
double  injectionX8coordinate = td.cloud().posP8_[0];
double  injectionX9coordinate = td.cloud().posP9_[0];
double  injectionX10coordinate = td.cloud().posP10_[0];
double  injectionY3coordinate = td.cloud().posP3_[1];
double  injectionY4coordinate = td.cloud().posP4_[1];
double  injectionY5coordinate = td.cloud().posP5_[1];
double  injectionY6coordinate = td.cloud().posP6_[1];
double  injectionY7coordinate = td.cloud().posP7_[1];
double  injectionY8coordinate = td.cloud().posP8_[1];
double  injectionY9coordinate = td.cloud().posP9_[1];
double  injectionY10coordinate = td.cloud().posP10_[1];
double  injectionX11coordinate = td.cloud().posP11_[0];
double  injectionX12coordinate = td.cloud().posP12_[0];
double  injectionX13coordinate = td.cloud().posP13_[0];
double  injectionX14coordinate = td.cloud().posP14_[0];
double  injectionX15coordinate = td.cloud().posP15_[0];
double  injectionX16coordinate = td.cloud().posP16_[0];
double  injectionX17coordinate = td.cloud().posP17_[0];
double  injectionX18coordinate = td.cloud().posP18_[0];
double  injectionX19coordinate = td.cloud().posP19_[0];
double  injectionX20coordinate = td.cloud().posP20_[0];
double  injectionX21coordinate = td.cloud().posP21_[0];
double  injectionX22coordinate = td.cloud().posP22_[0];
double  injectionX23coordinate = td.cloud().posP23_[0];
double  injectionX24coordinate = td.cloud().posP24_[0];
double  injectionX25coordinate = td.cloud().posP25_[0];
double  injectionX26coordinate = td.cloud().posP26_[0];
double  injectionX27coordinate = td.cloud().posP27_[0];
double  injectionX28coordinate = td.cloud().posP28_[0];
double  injectionX29coordinate = td.cloud().posP29_[0];
double  injectionX30coordinate = td.cloud().posP30_[0];
double  injectionX31coordinate = td.cloud().posP31_[0];
double  injectionX32coordinate = td.cloud().posP32_[0];
double  injectionX33coordinate = td.cloud().posP33_[0];
double  injectionX34coordinate = td.cloud().posP34_[0];
double  injectionX35coordinate = td.cloud().posP35_[0];
double  injectionX36coordinate = td.cloud().posP36_[0];
double  injectionX37coordinate = td.cloud().posP37_[0];
double  injectionX38coordinate = td.cloud().posP38_[0];
double  injectionX39coordinate = td.cloud().posP39_[0];
double  injectionX40coordinate = td.cloud().posP40_[0];
double  injectionX41coordinate = td.cloud().posP41_[0];
double  injectionX42coordinate = td.cloud().posP42_[0];
double  injectionX43coordinate = td.cloud().posP43_[0];
double  injectionX44coordinate = td.cloud().posP44_[0];
double  injectionX45coordinate = td.cloud().posP45_[0];
double  injectionX46coordinate = td.cloud().posP46_[0];
double  injectionX47coordinate = td.cloud().posP47_[0];
double  injectionX48coordinate = td.cloud().posP48_[0];
double  injectionX49coordinate = td.cloud().posP49_[0];
double  injectionX50coordinate = td.cloud().posP50_[0];
double  injectionX51coordinate = td.cloud().posP51_[0];
double  injectionX52coordinate = td.cloud().posP52_[0];
double  injectionX53coordinate = td.cloud().posP53_[0];
double  injectionX54coordinate = td.cloud().posP54_[0];
double  injectionX55coordinate = td.cloud().posP55_[0];
double  injectionX56coordinate = td.cloud().posP56_[0];
double  injectionX57coordinate = td.cloud().posP57_[0];
double  injectionX58coordinate = td.cloud().posP58_[0];
double  injectionX59coordinate = td.cloud().posP59_[0];
double  injectionX60coordinate = td.cloud().posP60_[0];
double  injectionX61coordinate = td.cloud().posP61_[0];
double  injectionX62coordinate = td.cloud().posP62_[0];
double  injectionX63coordinate = td.cloud().posP63_[0];
double  injectionX64coordinate = td.cloud().posP64_[0];
double  injectionX65coordinate = td.cloud().posP65_[0];
double  injectionX66coordinate = td.cloud().posP66_[0];
double  injectionX67coordinate = td.cloud().posP67_[0];
double  injectionX68coordinate = td.cloud().posP68_[0];
double  injectionX69coordinate = td.cloud().posP69_[0];
double  injectionX70coordinate = td.cloud().posP70_[0];
double  injectionX71coordinate = td.cloud().posP71_[0];
double  injectionX72coordinate = td.cloud().posP72_[0];
double  injectionX73coordinate = td.cloud().posP73_[0];
double  injectionX74coordinate = td.cloud().posP74_[0];
double  injectionX75coordinate = td.cloud().posP75_[0];
double  injectionX76coordinate = td.cloud().posP76_[0];
double  injectionX77coordinate = td.cloud().posP77_[0];
double  injectionX78coordinate = td.cloud().posP78_[0];
double  injectionX79coordinate = td.cloud().posP79_[0];
double  injectionX80coordinate = td.cloud().posP80_[0];
double  injectionX81coordinate = td.cloud().posP81_[0];
double  injectionX82coordinate = td.cloud().posP82_[0];
double  injectionX83coordinate = td.cloud().posP83_[0];
double  injectionX84coordinate = td.cloud().posP84_[0];
double  injectionX85coordinate = td.cloud().posP85_[0];
double  injectionX86coordinate = td.cloud().posP86_[0];
double  injectionX87coordinate = td.cloud().posP87_[0];
double  injectionX88coordinate = td.cloud().posP88_[0];
double  injectionX89coordinate = td.cloud().posP89_[0];
double  injectionX90coordinate = td.cloud().posP90_[0];
double  injectionX91coordinate = td.cloud().posP91_[0];
double  injectionX92coordinate = td.cloud().posP92_[0];
double  injectionX93coordinate = td.cloud().posP93_[0];
double  injectionX94coordinate = td.cloud().posP94_[0];
double  injectionX95coordinate = td.cloud().posP95_[0];
double  injectionX96coordinate = td.cloud().posP96_[0];
double  injectionX97coordinate = td.cloud().posP97_[0];
double  injectionX98coordinate = td.cloud().posP98_[0];
double  injectionX99coordinate = td.cloud().posP99_[0];
double  injectionX100coordinate = td.cloud().posP100_[0];


int  endP =  mesh_.points().size();
const point& firstPoint =mesh_.points()[0];
double xPartFirst =firstPoint[0];
const  point& lastPoint=mesh_.points()[999]; //need to add in by hand the point to take here, careful of this whenever the mesh resolution changes:including in the case of a serial run!!
double xPartLast =lastPoint[0];

Info<< "injection 1 x coord= " << injectionX1coordinate<< endl;
Info<< "injection 1 y coord= " << injectionY1coordinate<< endl;
Info<< "injection 2 x coord= " << injectionX2coordinate<< endl;
Info<< "injection 2 y coord= " << injectionY2coordinate<< endl;
Info<< "injection 3 x coord= " << injectionX3coordinate<< endl;
Info<< "injection 3 y coord= " << injectionY3coordinate<< endl;
Info<< "injection 4 x coord= " << injectionX4coordinate<< endl;
Info<< "injection 4 y coord= " << injectionY4coordinate<< endl;
Info<< "injection 5 x coord= " << injectionX5coordinate<< endl;
Info<< "injection 5 y coord= " << injectionY5coordinate<< endl;
Info<< "injection 6 x coord= " << injectionX6coordinate<< endl;
Info<< "injection 6 y coord= " << injectionY6coordinate<< endl;
Info<< "injection 7 x coord= " << injectionX7coordinate<< endl;
Info<< "injection 7 y coord= " << injectionY7coordinate<< endl;
Info<< "injection 8 x coord= " << injectionX8coordinate<< endl;
Info<< "injection 8 y coord= " << injectionY8coordinate<< endl;
Info<< "injection 9 x coord= " << injectionX9coordinate<< endl;
Info<< "injection 9 y coord= " << injectionY9coordinate<< endl;
Info<< "injection 10 x coord= " << injectionX10coordinate<< endl;
Info<< "injection 10 y coord= " << injectionY10coordinate<< endl;




Info<< "size= " << endP << endl;
Pout<< "firstPoint= " << xPartFirst<< endl;
Pout<< "lastPoint1= " << lastPoint<< endl;

  mesh_.findCellFacePt(td.cloud().posP1_, cellI, tetFaceI, tetPtI); 

if( (  injectionX1coordinate > xPartFirst ) && (  injectionX1coordinate <  xPartLast) )
{ Pout << "Particle 1 Injected on this processor" << endl;
  solidParticle* ptr1=new solidParticle(mesh_,td.cloud().posP1_,cellI,tetFaceI,tetPtI,td.cloud().dP1_,td.cloud().UP1_);
  Foam::Cloud<solidParticle>::addParticle(ptr1);
}

if((  injectionX2coordinate > xPartFirst ) && (  injectionX2coordinate <  xPartLast)) 
{ Pout << "Particle 2 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP2_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr2=new solidParticle(mesh_,td.cloud().posP2_,cellI,tetFaceI,tetPtI,td.cloud().dP2_,td.cloud().UP2_);
 Foam::Cloud<solidParticle>::addParticle(ptr2); }

if((  injectionX3coordinate > xPartFirst ) && (  injectionX3coordinate <  xPartLast)) 
{ Pout << "Particle 3 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP3_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr3=new solidParticle(mesh_,td.cloud().posP3_,cellI,tetFaceI,tetPtI,td.cloud().dP3_,td.cloud().UP3_);
 Foam::Cloud<solidParticle>::addParticle(ptr3); }

if((  injectionX4coordinate > xPartFirst ) && (  injectionX4coordinate <  xPartLast)) 
{ Pout << "Particle 4 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP4_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr4=new solidParticle(mesh_,td.cloud().posP4_,cellI,tetFaceI,tetPtI,td.cloud().dP4_,td.cloud().UP4_);
 Foam::Cloud<solidParticle>::addParticle(ptr4); }

if((  injectionX5coordinate > xPartFirst ) && (  injectionX5coordinate <  xPartLast)) 
{ Pout << "Particle 5 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP5_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr5=new solidParticle(mesh_,td.cloud().posP5_,cellI,tetFaceI,tetPtI,td.cloud().dP5_,td.cloud().UP5_);
 Foam::Cloud<solidParticle>::addParticle(ptr5); }

if((  injectionX6coordinate > xPartFirst ) && (  injectionX6coordinate <  xPartLast)) 
{ Pout << "Particle 6 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP6_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr6=new solidParticle(mesh_,td.cloud().posP6_,cellI,tetFaceI,tetPtI,td.cloud().dP6_,td.cloud().UP6_);
 Foam::Cloud<solidParticle>::addParticle(ptr6); }

if((  injectionX7coordinate > xPartFirst ) && (  injectionX7coordinate <  xPartLast)) 
{ Pout << "Particle 7 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP7_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr7=new solidParticle(mesh_,td.cloud().posP7_,cellI,tetFaceI,tetPtI,td.cloud().dP7_,td.cloud().UP7_);
 Foam::Cloud<solidParticle>::addParticle(ptr7); }

if((  injectionX8coordinate > xPartFirst ) && (  injectionX8coordinate <  xPartLast)) 
{ Pout << "Particle 8 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP8_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr8=new solidParticle(mesh_,td.cloud().posP8_,cellI,tetFaceI,tetPtI,td.cloud().dP8_,td.cloud().UP8_);
 Foam::Cloud<solidParticle>::addParticle(ptr8); }

if((  injectionX9coordinate > xPartFirst ) && (  injectionX9coordinate <  xPartLast)) 
{ Pout << "Particle 9 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP9_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr9=new solidParticle(mesh_,td.cloud().posP9_,cellI,tetFaceI,tetPtI,td.cloud().dP9_,td.cloud().UP9_);
 Foam::Cloud<solidParticle>::addParticle(ptr9); }

if((  injectionX10coordinate > xPartFirst ) && (  injectionX10coordinate <  xPartLast)) 
{ Pout << "Particle 10 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP10_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr10=new solidParticle(mesh_,td.cloud().posP10_,cellI,tetFaceI,tetPtI,td.cloud().dP10_,td.cloud().UP10_);
 Foam::Cloud<solidParticle>::addParticle(ptr10); }

if((  injectionX11coordinate > xPartFirst ) && (  injectionX11coordinate <  xPartLast)) 
{ Pout << "Particle 11 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP11_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr11=new solidParticle(mesh_,td.cloud().posP11_,cellI,tetFaceI,tetPtI,td.cloud().dP11_,td.cloud().UP11_);
 Foam::Cloud<solidParticle>::addParticle(ptr11); }

if((  injectionX12coordinate > xPartFirst ) && (  injectionX12coordinate <  xPartLast)) 
{ Pout << "Particle 12 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP12_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr12=new solidParticle(mesh_,td.cloud().posP12_,cellI,tetFaceI,tetPtI,td.cloud().dP12_,td.cloud().UP12_);
 Foam::Cloud<solidParticle>::addParticle(ptr12); }

if((  injectionX13coordinate > xPartFirst ) && (  injectionX13coordinate <  xPartLast)) 
{ Pout << "Particle 13 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP13_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr13=new solidParticle(mesh_,td.cloud().posP13_,cellI,tetFaceI,tetPtI,td.cloud().dP13_,td.cloud().UP13_);
 Foam::Cloud<solidParticle>::addParticle(ptr13); }

if((  injectionX14coordinate > xPartFirst ) && (  injectionX14coordinate <  xPartLast)) 
{ Pout << "Particle 14 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP14_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr14=new solidParticle(mesh_,td.cloud().posP14_,cellI,tetFaceI,tetPtI,td.cloud().dP14_,td.cloud().UP14_);
 Foam::Cloud<solidParticle>::addParticle(ptr14); }

if((  injectionX15coordinate > xPartFirst ) && (  injectionX15coordinate <  xPartLast)) 
{ Pout << "Particle 15 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP15_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr15=new solidParticle(mesh_,td.cloud().posP15_,cellI,tetFaceI,tetPtI,td.cloud().dP15_,td.cloud().UP15_);
 Foam::Cloud<solidParticle>::addParticle(ptr15); }

if((  injectionX16coordinate > xPartFirst ) && (  injectionX16coordinate <  xPartLast)) 
{ Pout << "Particle 16 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP16_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr16=new solidParticle(mesh_,td.cloud().posP16_,cellI,tetFaceI,tetPtI,td.cloud().dP16_,td.cloud().UP16_);
 Foam::Cloud<solidParticle>::addParticle(ptr16); 
Pout << "All Good" << endl;}

if((  injectionX17coordinate > xPartFirst ) && (  injectionX17coordinate <  xPartLast)) 
{ Pout << "Particle 17 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP17_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr17=new solidParticle(mesh_,td.cloud().posP17_,cellI,tetFaceI,tetPtI,td.cloud().dP17_,td.cloud().UP17_);
 Foam::Cloud<solidParticle>::addParticle(ptr17); }

if((  injectionX18coordinate > xPartFirst ) && (  injectionX18coordinate <  xPartLast)) 
{ Pout << "Particle 18 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP18_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr18=new solidParticle(mesh_,td.cloud().posP18_,cellI,tetFaceI,tetPtI,td.cloud().dP18_,td.cloud().UP18_);
 Foam::Cloud<solidParticle>::addParticle(ptr18); }

if((  injectionX19coordinate > xPartFirst ) && (  injectionX19coordinate <  xPartLast)) 
{ Pout << "Particle 19 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP19_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr19=new solidParticle(mesh_,td.cloud().posP19_,cellI,tetFaceI,tetPtI,td.cloud().dP19_,td.cloud().UP19_);
 Foam::Cloud<solidParticle>::addParticle(ptr19); }

if((  injectionX20coordinate > xPartFirst ) && (  injectionX20coordinate <  xPartLast)) 
{ Pout << "Particle 20 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP20_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr20=new solidParticle(mesh_,td.cloud().posP20_,cellI,tetFaceI,tetPtI,td.cloud().dP20_,td.cloud().UP20_);
 Foam::Cloud<solidParticle>::addParticle(ptr20); }

if((  injectionX21coordinate > xPartFirst ) && (  injectionX21coordinate <  xPartLast)) 
{ Pout << "Particle 21 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP21_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr21=new solidParticle(mesh_,td.cloud().posP21_,cellI,tetFaceI,tetPtI,td.cloud().dP21_,td.cloud().UP21_);
 Foam::Cloud<solidParticle>::addParticle(ptr21); }

if((  injectionX22coordinate > xPartFirst ) && (  injectionX22coordinate <  xPartLast)) 
{ Pout << "Particle 22 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP22_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr22=new solidParticle(mesh_,td.cloud().posP22_,cellI,tetFaceI,tetPtI,td.cloud().dP22_,td.cloud().UP22_);
 Foam::Cloud<solidParticle>::addParticle(ptr22); }

if((  injectionX23coordinate > xPartFirst ) && (  injectionX23coordinate <  xPartLast)) 
{ Pout << "Particle 23 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP23_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr23=new solidParticle(mesh_,td.cloud().posP23_,cellI,tetFaceI,tetPtI,td.cloud().dP23_,td.cloud().UP23_);
 Foam::Cloud<solidParticle>::addParticle(ptr23); }

if((  injectionX24coordinate > xPartFirst ) && (  injectionX24coordinate <  xPartLast)) 
{ Pout << "Particle 24 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP24_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr24=new solidParticle(mesh_,td.cloud().posP24_,cellI,tetFaceI,tetPtI,td.cloud().dP24_,td.cloud().UP24_);
 Foam::Cloud<solidParticle>::addParticle(ptr24); }

if((  injectionX25coordinate > xPartFirst ) && (  injectionX25coordinate <  xPartLast)) 
{ Pout << "Particle 25 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP25_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr25=new solidParticle(mesh_,td.cloud().posP25_,cellI,tetFaceI,tetPtI,td.cloud().dP25_,td.cloud().UP25_);
 Foam::Cloud<solidParticle>::addParticle(ptr25); }

if((  injectionX26coordinate > xPartFirst ) && (  injectionX26coordinate <  xPartLast)) 
{ Pout << "Particle 26 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP26_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr26=new solidParticle(mesh_,td.cloud().posP26_,cellI,tetFaceI,tetPtI,td.cloud().dP26_,td.cloud().UP26_);
 Foam::Cloud<solidParticle>::addParticle(ptr26); }

if((  injectionX27coordinate > xPartFirst ) && (  injectionX27coordinate <  xPartLast)) 
{ Pout << "Particle 27 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP27_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr27=new solidParticle(mesh_,td.cloud().posP27_,cellI,tetFaceI,tetPtI,td.cloud().dP27_,td.cloud().UP27_);
 Foam::Cloud<solidParticle>::addParticle(ptr27); }

if((  injectionX28coordinate > xPartFirst ) && (  injectionX28coordinate <  xPartLast)) 
{ Pout << "Particle 28 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP28_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr28=new solidParticle(mesh_,td.cloud().posP28_,cellI,tetFaceI,tetPtI,td.cloud().dP28_,td.cloud().UP28_);
 Foam::Cloud<solidParticle>::addParticle(ptr28); }

if((  injectionX29coordinate > xPartFirst ) && (  injectionX29coordinate <  xPartLast)) 
{ Pout << "Particle 29 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP29_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr29=new solidParticle(mesh_,td.cloud().posP29_,cellI,tetFaceI,tetPtI,td.cloud().dP29_,td.cloud().UP29_);
 Foam::Cloud<solidParticle>::addParticle(ptr29); }

if((  injectionX30coordinate > xPartFirst ) && (  injectionX30coordinate <  xPartLast)) 
{ Pout << "Particle 30 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP30_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr30=new solidParticle(mesh_,td.cloud().posP30_,cellI,tetFaceI,tetPtI,td.cloud().dP30_,td.cloud().UP30_);
 Foam::Cloud<solidParticle>::addParticle(ptr30); }

if((  injectionX31coordinate > xPartFirst ) && (  injectionX31coordinate <  xPartLast)) 
{ Pout << "Particle 31 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP31_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr31=new solidParticle(mesh_,td.cloud().posP31_,cellI,tetFaceI,tetPtI,td.cloud().dP31_,td.cloud().UP31_);
 Foam::Cloud<solidParticle>::addParticle(ptr31); }

if((  injectionX32coordinate > xPartFirst ) && (  injectionX32coordinate <  xPartLast)) 
{ Pout << "Particle 32 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP32_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr32=new solidParticle(mesh_,td.cloud().posP32_,cellI,tetFaceI,tetPtI,td.cloud().dP32_,td.cloud().UP32_);
 Foam::Cloud<solidParticle>::addParticle(ptr32); }

if((  injectionX33coordinate > xPartFirst ) && (  injectionX33coordinate <  xPartLast)) 
{ Pout << "Particle 33 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP33_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr33=new solidParticle(mesh_,td.cloud().posP33_,cellI,tetFaceI,tetPtI,td.cloud().dP33_,td.cloud().UP33_);
 Foam::Cloud<solidParticle>::addParticle(ptr33); }

if((  injectionX34coordinate > xPartFirst ) && (  injectionX34coordinate <  xPartLast)) 
{ Pout << "Particle 34 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP34_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr34=new solidParticle(mesh_,td.cloud().posP34_,cellI,tetFaceI,tetPtI,td.cloud().dP34_,td.cloud().UP34_);
 Foam::Cloud<solidParticle>::addParticle(ptr34); }

if((  injectionX35coordinate > xPartFirst ) && (  injectionX35coordinate <  xPartLast)) 
{ Pout << "Particle 35 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP35_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr35=new solidParticle(mesh_,td.cloud().posP35_,cellI,tetFaceI,tetPtI,td.cloud().dP35_,td.cloud().UP35_);
 Foam::Cloud<solidParticle>::addParticle(ptr35); }

if((  injectionX36coordinate > xPartFirst ) && (  injectionX36coordinate <  xPartLast)) 
{ Pout << "Particle 36 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP36_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr36=new solidParticle(mesh_,td.cloud().posP36_,cellI,tetFaceI,tetPtI,td.cloud().dP36_,td.cloud().UP36_);
 Foam::Cloud<solidParticle>::addParticle(ptr36); }

if((  injectionX37coordinate > xPartFirst ) && (  injectionX37coordinate <  xPartLast)) 
{ Pout << "Particle 37 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP37_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr37=new solidParticle(mesh_,td.cloud().posP37_,cellI,tetFaceI,tetPtI,td.cloud().dP37_,td.cloud().UP37_);
 Foam::Cloud<solidParticle>::addParticle(ptr37); }

if((  injectionX38coordinate > xPartFirst ) && (  injectionX38coordinate <  xPartLast)) 
{ Pout << "Particle 38 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP38_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr38=new solidParticle(mesh_,td.cloud().posP38_,cellI,tetFaceI,tetPtI,td.cloud().dP38_,td.cloud().UP38_);
 Foam::Cloud<solidParticle>::addParticle(ptr38); }

if((  injectionX39coordinate > xPartFirst ) && (  injectionX39coordinate <  xPartLast)) 
{ Pout << "Particle 39 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP39_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr39=new solidParticle(mesh_,td.cloud().posP39_,cellI,tetFaceI,tetPtI,td.cloud().dP39_,td.cloud().UP39_);
 Foam::Cloud<solidParticle>::addParticle(ptr39); }

if((  injectionX40coordinate > xPartFirst ) && (  injectionX40coordinate <  xPartLast)) 
{ Pout << "Particle 40 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP40_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr40=new solidParticle(mesh_,td.cloud().posP40_,cellI,tetFaceI,tetPtI,td.cloud().dP40_,td.cloud().UP40_);
 Foam::Cloud<solidParticle>::addParticle(ptr40); }

if((  injectionX41coordinate > xPartFirst ) && (  injectionX41coordinate <  xPartLast)) 
{ Pout << "Particle 41 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP41_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr41=new solidParticle(mesh_,td.cloud().posP41_,cellI,tetFaceI,tetPtI,td.cloud().dP41_,td.cloud().UP41_);
 Foam::Cloud<solidParticle>::addParticle(ptr41); }

if((  injectionX42coordinate > xPartFirst ) && (  injectionX42coordinate <  xPartLast)) 
{ Pout << "Particle 42 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP42_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr42=new solidParticle(mesh_,td.cloud().posP42_,cellI,tetFaceI,tetPtI,td.cloud().dP42_,td.cloud().UP42_);
 Foam::Cloud<solidParticle>::addParticle(ptr42); }

if((  injectionX43coordinate > xPartFirst ) && (  injectionX43coordinate <  xPartLast)) 
{ Pout << "Particle 43 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP43_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr43=new solidParticle(mesh_,td.cloud().posP43_,cellI,tetFaceI,tetPtI,td.cloud().dP43_,td.cloud().UP43_);
 Foam::Cloud<solidParticle>::addParticle(ptr43); }

if((  injectionX44coordinate > xPartFirst ) && (  injectionX44coordinate <  xPartLast)) 
{ Pout << "Particle 44 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP44_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr44=new solidParticle(mesh_,td.cloud().posP44_,cellI,tetFaceI,tetPtI,td.cloud().dP44_,td.cloud().UP44_);
 Foam::Cloud<solidParticle>::addParticle(ptr44); }

if((  injectionX45coordinate > xPartFirst ) && (  injectionX45coordinate <  xPartLast)) 
{ Pout << "Particle 45 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP45_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr45=new solidParticle(mesh_,td.cloud().posP45_,cellI,tetFaceI,tetPtI,td.cloud().dP45_,td.cloud().UP45_);
 Foam::Cloud<solidParticle>::addParticle(ptr45); }

if((  injectionX46coordinate > xPartFirst ) && (  injectionX46coordinate <  xPartLast)) 
{ Pout << "Particle 46 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP46_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr46=new solidParticle(mesh_,td.cloud().posP46_,cellI,tetFaceI,tetPtI,td.cloud().dP46_,td.cloud().UP46_);
 Foam::Cloud<solidParticle>::addParticle(ptr46); }

if((  injectionX47coordinate > xPartFirst ) && (  injectionX47coordinate <  xPartLast)) 
{ Pout << "Particle 47 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP47_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr47=new solidParticle(mesh_,td.cloud().posP47_,cellI,tetFaceI,tetPtI,td.cloud().dP47_,td.cloud().UP47_);
 Foam::Cloud<solidParticle>::addParticle(ptr47); }

if((  injectionX48coordinate > xPartFirst ) && (  injectionX48coordinate <  xPartLast)) 
{ Pout << "Particle 48 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP48_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr48=new solidParticle(mesh_,td.cloud().posP48_,cellI,tetFaceI,tetPtI,td.cloud().dP48_,td.cloud().UP48_);
 Foam::Cloud<solidParticle>::addParticle(ptr48); }

if((  injectionX49coordinate > xPartFirst ) && (  injectionX49coordinate <  xPartLast)) 
{ Pout << "Particle 49 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP49_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr49=new solidParticle(mesh_,td.cloud().posP49_,cellI,tetFaceI,tetPtI,td.cloud().dP49_,td.cloud().UP49_);
 Foam::Cloud<solidParticle>::addParticle(ptr49); }

if((  injectionX50coordinate > xPartFirst ) && (  injectionX50coordinate <  xPartLast)) 
{ Pout << "Particle 50 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP50_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr50=new solidParticle(mesh_,td.cloud().posP50_,cellI,tetFaceI,tetPtI,td.cloud().dP50_,td.cloud().UP50_);
 Foam::Cloud<solidParticle>::addParticle(ptr50); }





if((  injectionX51coordinate > xPartFirst ) && (  injectionX51coordinate <  xPartLast)) 
{ Pout << "Particle 51 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP51_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr51=new solidParticle(mesh_,td.cloud().posP51_,cellI,tetFaceI,tetPtI,td.cloud().dP51_,td.cloud().UP51_);
 Foam::Cloud<solidParticle>::addParticle(ptr51); }

if((  injectionX52coordinate > xPartFirst ) && (  injectionX52coordinate <  xPartLast)) 
{ Pout << "Particle 52 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP52_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr52=new solidParticle(mesh_,td.cloud().posP52_,cellI,tetFaceI,tetPtI,td.cloud().dP52_,td.cloud().UP52_);
 Foam::Cloud<solidParticle>::addParticle(ptr52); }

if((  injectionX53coordinate > xPartFirst ) && (  injectionX53coordinate <  xPartLast)) 
{ Pout << "Particle 53 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP53_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr53=new solidParticle(mesh_,td.cloud().posP53_,cellI,tetFaceI,tetPtI,td.cloud().dP53_,td.cloud().UP53_);
 Foam::Cloud<solidParticle>::addParticle(ptr53); }

if((  injectionX54coordinate > xPartFirst ) && (  injectionX54coordinate <  xPartLast)) 
{ Pout << "Particle 54 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP54_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr54=new solidParticle(mesh_,td.cloud().posP54_,cellI,tetFaceI,tetPtI,td.cloud().dP54_,td.cloud().UP54_);
 Foam::Cloud<solidParticle>::addParticle(ptr54); }

if((  injectionX55coordinate > xPartFirst ) && (  injectionX55coordinate <  xPartLast)) 
{ Pout << "Particle 55 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP55_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr55=new solidParticle(mesh_,td.cloud().posP55_,cellI,tetFaceI,tetPtI,td.cloud().dP55_,td.cloud().UP55_);
 Foam::Cloud<solidParticle>::addParticle(ptr55); }

if((  injectionX56coordinate > xPartFirst ) && (  injectionX56coordinate <  xPartLast)) 
{ Pout << "Particle 56 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP56_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr56=new solidParticle(mesh_,td.cloud().posP56_,cellI,tetFaceI,tetPtI,td.cloud().dP56_,td.cloud().UP56_);
 Foam::Cloud<solidParticle>::addParticle(ptr56); }

if((  injectionX57coordinate > xPartFirst ) && (  injectionX57coordinate <  xPartLast)) 
{ Pout << "Particle 57 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP57_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr57=new solidParticle(mesh_,td.cloud().posP57_,cellI,tetFaceI,tetPtI,td.cloud().dP57_,td.cloud().UP57_);
 Foam::Cloud<solidParticle>::addParticle(ptr57); }

if((  injectionX58coordinate > xPartFirst ) && (  injectionX58coordinate <  xPartLast)) 
{ Pout << "Particle 58 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP58_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr58=new solidParticle(mesh_,td.cloud().posP58_,cellI,tetFaceI,tetPtI,td.cloud().dP58_,td.cloud().UP58_);
 Foam::Cloud<solidParticle>::addParticle(ptr58); }

if((  injectionX59coordinate > xPartFirst ) && (  injectionX59coordinate <  xPartLast)) 
{ Pout << "Particle 59 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP59_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr59=new solidParticle(mesh_,td.cloud().posP59_,cellI,tetFaceI,tetPtI,td.cloud().dP59_,td.cloud().UP59_);
 Foam::Cloud<solidParticle>::addParticle(ptr59); }

if((  injectionX60coordinate > xPartFirst ) && (  injectionX60coordinate <  xPartLast)) 
{ Pout << "Particle 60 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP60_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr60=new solidParticle(mesh_,td.cloud().posP60_,cellI,tetFaceI,tetPtI,td.cloud().dP60_,td.cloud().UP60_);
 Foam::Cloud<solidParticle>::addParticle(ptr60); }

if((  injectionX61coordinate > xPartFirst ) && (  injectionX61coordinate <  xPartLast)) 
{ Pout << "Particle 61 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP61_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr61=new solidParticle(mesh_,td.cloud().posP61_,cellI,tetFaceI,tetPtI,td.cloud().dP61_,td.cloud().UP61_);
 Foam::Cloud<solidParticle>::addParticle(ptr61); }

if((  injectionX62coordinate > xPartFirst ) && (  injectionX62coordinate <  xPartLast)) 
{ Pout << "Particle 62 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP62_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr62=new solidParticle(mesh_,td.cloud().posP62_,cellI,tetFaceI,tetPtI,td.cloud().dP62_,td.cloud().UP62_);
 Foam::Cloud<solidParticle>::addParticle(ptr62); }

if((  injectionX63coordinate > xPartFirst ) && (  injectionX63coordinate <  xPartLast)) 
{ Pout << "Particle 63 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP63_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr63=new solidParticle(mesh_,td.cloud().posP63_,cellI,tetFaceI,tetPtI,td.cloud().dP63_,td.cloud().UP63_);
 Foam::Cloud<solidParticle>::addParticle(ptr63); }

if((  injectionX64coordinate > xPartFirst ) && (  injectionX64coordinate <  xPartLast)) 
{ Pout << "Particle 64 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP64_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr64=new solidParticle(mesh_,td.cloud().posP64_,cellI,tetFaceI,tetPtI,td.cloud().dP64_,td.cloud().UP64_);
 Foam::Cloud<solidParticle>::addParticle(ptr64); }

if((  injectionX65coordinate > xPartFirst ) && (  injectionX65coordinate <  xPartLast)) 
{ Pout << "Particle 65 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP65_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr65=new solidParticle(mesh_,td.cloud().posP65_,cellI,tetFaceI,tetPtI,td.cloud().dP65_,td.cloud().UP65_);
 Foam::Cloud<solidParticle>::addParticle(ptr65); }

if((  injectionX66coordinate > xPartFirst ) && (  injectionX66coordinate <  xPartLast)) 
{ Pout << "Particle 66 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP66_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr66=new solidParticle(mesh_,td.cloud().posP66_,cellI,tetFaceI,tetPtI,td.cloud().dP66_,td.cloud().UP66_);
 Foam::Cloud<solidParticle>::addParticle(ptr66); }

if((  injectionX67coordinate > xPartFirst ) && (  injectionX67coordinate <  xPartLast)) 
{ Pout << "Particle 67 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP67_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr67=new solidParticle(mesh_,td.cloud().posP67_,cellI,tetFaceI,tetPtI,td.cloud().dP67_,td.cloud().UP67_);
 Foam::Cloud<solidParticle>::addParticle(ptr67); }

if((  injectionX68coordinate > xPartFirst ) && (  injectionX68coordinate <  xPartLast)) 
{ Pout << "Particle 68 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP68_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr68=new solidParticle(mesh_,td.cloud().posP68_,cellI,tetFaceI,tetPtI,td.cloud().dP68_,td.cloud().UP68_);
 Foam::Cloud<solidParticle>::addParticle(ptr68); }

if((  injectionX69coordinate > xPartFirst ) && (  injectionX69coordinate <  xPartLast)) 
{ Pout << "Particle 69 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP69_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr69=new solidParticle(mesh_,td.cloud().posP69_,cellI,tetFaceI,tetPtI,td.cloud().dP69_,td.cloud().UP69_);
 Foam::Cloud<solidParticle>::addParticle(ptr69); }

if((  injectionX70coordinate > xPartFirst ) && (  injectionX70coordinate <  xPartLast)) 
{ Pout << "Particle 70 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP70_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr70=new solidParticle(mesh_,td.cloud().posP70_,cellI,tetFaceI,tetPtI,td.cloud().dP70_,td.cloud().UP70_);
 Foam::Cloud<solidParticle>::addParticle(ptr70); }

if((  injectionX71coordinate > xPartFirst ) && (  injectionX71coordinate <  xPartLast)) 
{ Pout << "Particle 71 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP71_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr71=new solidParticle(mesh_,td.cloud().posP71_,cellI,tetFaceI,tetPtI,td.cloud().dP71_,td.cloud().UP71_);
 Foam::Cloud<solidParticle>::addParticle(ptr71); }

if((  injectionX72coordinate > xPartFirst ) && (  injectionX72coordinate <  xPartLast)) 
{ Pout << "Particle 72 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP72_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr72=new solidParticle(mesh_,td.cloud().posP72_,cellI,tetFaceI,tetPtI,td.cloud().dP72_,td.cloud().UP72_);
 Foam::Cloud<solidParticle>::addParticle(ptr72); }

if((  injectionX73coordinate > xPartFirst ) && (  injectionX73coordinate <  xPartLast)) 
{ Pout << "Particle 73 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP73_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr73=new solidParticle(mesh_,td.cloud().posP73_,cellI,tetFaceI,tetPtI,td.cloud().dP73_,td.cloud().UP73_);
 Foam::Cloud<solidParticle>::addParticle(ptr73); }

if((  injectionX74coordinate > xPartFirst ) && (  injectionX74coordinate <  xPartLast)) 
{ Pout << "Particle 74 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP74_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr74=new solidParticle(mesh_,td.cloud().posP74_,cellI,tetFaceI,tetPtI,td.cloud().dP74_,td.cloud().UP74_);
 Foam::Cloud<solidParticle>::addParticle(ptr74); }

if((  injectionX75coordinate > xPartFirst ) && (  injectionX75coordinate <  xPartLast)) 
{ Pout << "Particle 75 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP75_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr75=new solidParticle(mesh_,td.cloud().posP75_,cellI,tetFaceI,tetPtI,td.cloud().dP75_,td.cloud().UP75_);
 Foam::Cloud<solidParticle>::addParticle(ptr75); }

if((  injectionX76coordinate > xPartFirst ) && (  injectionX76coordinate <  xPartLast)) 
{ Pout << "Particle 76 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP76_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr76=new solidParticle(mesh_,td.cloud().posP76_,cellI,tetFaceI,tetPtI,td.cloud().dP76_,td.cloud().UP76_);
 Foam::Cloud<solidParticle>::addParticle(ptr76); }

if((  injectionX77coordinate > xPartFirst ) && (  injectionX77coordinate <  xPartLast)) 
{ Pout << "Particle 77 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP77_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr77=new solidParticle(mesh_,td.cloud().posP77_,cellI,tetFaceI,tetPtI,td.cloud().dP77_,td.cloud().UP77_);
 Foam::Cloud<solidParticle>::addParticle(ptr77); }

if((  injectionX78coordinate > xPartFirst ) && (  injectionX78coordinate <  xPartLast)) 
{ Pout << "Particle 78 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP78_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr78=new solidParticle(mesh_,td.cloud().posP78_,cellI,tetFaceI,tetPtI,td.cloud().dP78_,td.cloud().UP78_);
 Foam::Cloud<solidParticle>::addParticle(ptr78); }

if((  injectionX79coordinate > xPartFirst ) && (  injectionX79coordinate <  xPartLast)) 
{ Pout << "Particle 79 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP79_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr79=new solidParticle(mesh_,td.cloud().posP79_,cellI,tetFaceI,tetPtI,td.cloud().dP79_,td.cloud().UP79_);
 Foam::Cloud<solidParticle>::addParticle(ptr79); }

if((  injectionX80coordinate > xPartFirst ) && (  injectionX80coordinate <  xPartLast)) 
{ Pout << "Particle 80 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP80_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr80=new solidParticle(mesh_,td.cloud().posP80_,cellI,tetFaceI,tetPtI,td.cloud().dP80_,td.cloud().UP80_);
 Foam::Cloud<solidParticle>::addParticle(ptr80); }

if((  injectionX81coordinate > xPartFirst ) && (  injectionX81coordinate <  xPartLast)) 
{ Pout << "Particle 81 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP81_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr81=new solidParticle(mesh_,td.cloud().posP81_,cellI,tetFaceI,tetPtI,td.cloud().dP81_,td.cloud().UP81_);
 Foam::Cloud<solidParticle>::addParticle(ptr81); }

if((  injectionX82coordinate > xPartFirst ) && (  injectionX82coordinate <  xPartLast)) 
{ Pout << "Particle 82 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP82_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr82=new solidParticle(mesh_,td.cloud().posP82_,cellI,tetFaceI,tetPtI,td.cloud().dP82_,td.cloud().UP82_);
 Foam::Cloud<solidParticle>::addParticle(ptr82); }

if((  injectionX83coordinate > xPartFirst ) && (  injectionX83coordinate <  xPartLast)) 
{ Pout << "Particle 83 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP83_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr83=new solidParticle(mesh_,td.cloud().posP83_,cellI,tetFaceI,tetPtI,td.cloud().dP83_,td.cloud().UP83_);
 Foam::Cloud<solidParticle>::addParticle(ptr83); }

if((  injectionX84coordinate > xPartFirst ) && (  injectionX84coordinate <  xPartLast)) 
{ Pout << "Particle 84 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP84_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr84=new solidParticle(mesh_,td.cloud().posP84_,cellI,tetFaceI,tetPtI,td.cloud().dP84_,td.cloud().UP84_);
 Foam::Cloud<solidParticle>::addParticle(ptr84); }

if((  injectionX85coordinate > xPartFirst ) && (  injectionX85coordinate <  xPartLast)) 
{ Pout << "Particle 85 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP85_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr85=new solidParticle(mesh_,td.cloud().posP85_,cellI,tetFaceI,tetPtI,td.cloud().dP85_,td.cloud().UP85_);
 Foam::Cloud<solidParticle>::addParticle(ptr85); }

if((  injectionX86coordinate > xPartFirst ) && (  injectionX86coordinate <  xPartLast)) 
{ Pout << "Particle 86 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP86_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr86=new solidParticle(mesh_,td.cloud().posP86_,cellI,tetFaceI,tetPtI,td.cloud().dP86_,td.cloud().UP86_);
 Foam::Cloud<solidParticle>::addParticle(ptr86); }

if((  injectionX87coordinate > xPartFirst ) && (  injectionX87coordinate <  xPartLast)) 
{ Pout << "Particle 87 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP87_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr87=new solidParticle(mesh_,td.cloud().posP87_,cellI,tetFaceI,tetPtI,td.cloud().dP87_,td.cloud().UP87_);
 Foam::Cloud<solidParticle>::addParticle(ptr87); }

if((  injectionX88coordinate > xPartFirst ) && (  injectionX88coordinate <  xPartLast)) 
{ Pout << "Particle 88 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP88_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr88=new solidParticle(mesh_,td.cloud().posP88_,cellI,tetFaceI,tetPtI,td.cloud().dP88_,td.cloud().UP88_);
 Foam::Cloud<solidParticle>::addParticle(ptr88); }

if((  injectionX89coordinate > xPartFirst ) && (  injectionX89coordinate <  xPartLast)) 
{ Pout << "Particle 89 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP89_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr89=new solidParticle(mesh_,td.cloud().posP89_,cellI,tetFaceI,tetPtI,td.cloud().dP89_,td.cloud().UP89_);
 Foam::Cloud<solidParticle>::addParticle(ptr89); }

if((  injectionX90coordinate > xPartFirst ) && (  injectionX90coordinate <  xPartLast)) 
{ Pout << "Particle 90 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP90_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr90=new solidParticle(mesh_,td.cloud().posP90_,cellI,tetFaceI,tetPtI,td.cloud().dP90_,td.cloud().UP90_);
 Foam::Cloud<solidParticle>::addParticle(ptr90); }

if((  injectionX91coordinate > xPartFirst ) && (  injectionX91coordinate <  xPartLast)) 
{ Pout << "Particle 91 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP91_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr91=new solidParticle(mesh_,td.cloud().posP91_,cellI,tetFaceI,tetPtI,td.cloud().dP91_,td.cloud().UP91_);
 Foam::Cloud<solidParticle>::addParticle(ptr91); }

if((  injectionX92coordinate > xPartFirst ) && (  injectionX92coordinate <  xPartLast)) 
{ Pout << "Particle 92 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP92_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr92=new solidParticle(mesh_,td.cloud().posP92_,cellI,tetFaceI,tetPtI,td.cloud().dP92_,td.cloud().UP92_);
 Foam::Cloud<solidParticle>::addParticle(ptr92); }

if((  injectionX93coordinate > xPartFirst ) && (  injectionX93coordinate <  xPartLast)) 
{ Pout << "Particle 93 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP93_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr93=new solidParticle(mesh_,td.cloud().posP93_,cellI,tetFaceI,tetPtI,td.cloud().dP93_,td.cloud().UP93_);
 Foam::Cloud<solidParticle>::addParticle(ptr93); }

if((  injectionX94coordinate > xPartFirst ) && (  injectionX94coordinate <  xPartLast)) 
{ Pout << "Particle 94 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP94_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr94=new solidParticle(mesh_,td.cloud().posP94_,cellI,tetFaceI,tetPtI,td.cloud().dP94_,td.cloud().UP94_);
 Foam::Cloud<solidParticle>::addParticle(ptr94); }

if((  injectionX95coordinate > xPartFirst ) && (  injectionX95coordinate <  xPartLast)) 
{ Pout << "Particle 95 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP95_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr95=new solidParticle(mesh_,td.cloud().posP95_,cellI,tetFaceI,tetPtI,td.cloud().dP95_,td.cloud().UP95_);
 Foam::Cloud<solidParticle>::addParticle(ptr95); }

if((  injectionX96coordinate > xPartFirst ) && (  injectionX96coordinate <  xPartLast)) 
{ Pout << "Particle 96 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP96_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr96=new solidParticle(mesh_,td.cloud().posP96_,cellI,tetFaceI,tetPtI,td.cloud().dP96_,td.cloud().UP96_);
 Foam::Cloud<solidParticle>::addParticle(ptr96); }

if((  injectionX97coordinate > xPartFirst ) && (  injectionX97coordinate <  xPartLast)) 
{ Pout << "Particle 97 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP97_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr97=new solidParticle(mesh_,td.cloud().posP97_,cellI,tetFaceI,tetPtI,td.cloud().dP97_,td.cloud().UP97_);
 Foam::Cloud<solidParticle>::addParticle(ptr97); }

if((  injectionX98coordinate > xPartFirst ) && (  injectionX98coordinate <  xPartLast)) 
{ Pout << "Particle 98 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP98_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr98=new solidParticle(mesh_,td.cloud().posP98_,cellI,tetFaceI,tetPtI,td.cloud().dP98_,td.cloud().UP98_);
 Foam::Cloud<solidParticle>::addParticle(ptr98); }

if((  injectionX99coordinate > xPartFirst ) && (  injectionX99coordinate <  xPartLast)) 
{ Pout << "Particle 99 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP99_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr99=new solidParticle(mesh_,td.cloud().posP99_,cellI,tetFaceI,tetPtI,td.cloud().dP99_,td.cloud().UP99_);
 Foam::Cloud<solidParticle>::addParticle(ptr99); }

if((  injectionX100coordinate > xPartFirst ) && (  injectionX100coordinate <  xPartLast)) 
{ Pout << "Particle 100 Injected on this processor" << endl;
  mesh_.findCellFacePt(td.cloud().posP100_, cellI, tetFaceI, tetPtI);
   solidParticle* ptr100=new solidParticle(mesh_,td.cloud().posP100_,cellI,tetFaceI,tetPtI,td.cloud().dP100_,td.cloud().UP100_);
 Foam::Cloud<solidParticle>::addParticle(ptr100); }


}
bool Foam::solidParticleCloud::hasWallImpactDistance() const
{
    return true;
}
// ************************************************************************* //
