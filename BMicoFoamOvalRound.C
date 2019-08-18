/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/


#include "fvMesh.H"
#include "motionSolver.H"
#include "velocityMotionSolver.H"
#include "primitivePatchInterpolation.H"
#include "fvCFD.H"
#include "pisoControl.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"

#include "CorrectPhi.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"



#include "fvMesh.H"
#include "motionSolver.H"
#include "velocityMotionSolver.H"
#include "primitivePatchInterpolation.H"

#include "solidParticle.C"
#include "solidParticleIO.C"
#include "solidParticleCloud.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
	#include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

	#include "solidParticle.H"
	#include "solidParticleCloud.H"



    pisoControl piso(mesh);

	solidParticleCloud particles(mesh);
	Info<< "Cloud size= "<< particles.size() <<endl;



    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "createControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

	#include "findBM.H"
	#include "findOval.H"
  #include "findRound.H"

	int intPhase=0;


	Info << "\nBM has been found. Begin read in of the coordinate matrices\n" << endl;

	std::fstream dimFile("dimensions_y.txt", std::ios_base::in);
	std::fstream coordsFileY("matrix_y.txt", std::ios_base::in);
	std::fstream coordsFileX("matrix_x.txt", std::ios_base::in);

	std::fstream dimFile2("dimensions_Oval.txt", std::ios_base::in);
	std::fstream ovalYCoords("ovalYcoords.txt", std::ios_base::in);
	std::fstream ovalXCoords("ovalXcoords.txt", std::ios_base::in);

	std::fstream dimFile3("dimensions_Round.txt", std::ios_base::in);
	std::fstream roundYCoords("roundYcoords.txt", std::ios_base::in);
	std::fstream roundXCoords("roundXcoords.txt", std::ios_base::in);
  	float a, b, c, d;

	dimFile >> a >> b;
	dimFile2 >> c >> d;
//dimensions of round and oval windows always the same

Info << "a= " << a << " (# of timesteps)" << endl;
Info << "b= " << b << " (# of BM grid points)" << endl;
Info << "c= " << c << " (# of timesteps)" << endl;
Info << "d= " << d << " (# of Oval Window grid points)" << endl;
	std::vector<std::vector<double> > matrixX(a, std::vector<double>(b));
	std::vector<std::vector<double> > matrixY(a, std::vector<double>(b));
	std::vector<std::vector<double> > matrixOvalX(c, std::vector<double>(d));
	std::vector<std::vector<double> > matrixOvalY(c, std::vector<double>(d));
  std::vector<std::vector<double> > matrixRoundX(c, std::vector<double>(d));
	std::vector<std::vector<double> > matrixRoundY(c, std::vector<double>(d));

if(coordsFileX.is_open())
{

	Info << "File matrix_x opened! Reading data from file into coordsFileX" << endl;
for (int x1 = 0; x1 < a; x1++)
	{

    for (int y1 = 0; y1 < b; y1++)
	{

        coordsFileX >> matrixX[x1][y1];

    }
  }
				}

else
	{
		Info << "File matrix_x not found! Abort" << endl;
		exit(1);
	}

if(coordsFileY.is_open())
{

	//Info << "File " << coordsFileY.c_str() <<" opened! Reading data from file into matrix_y" << endl;
for (int x2 = 0; x2 < a; x2++)
	{

    for (int y2 = 0; y2 < b; y2++)
	{

        coordsFileY >> matrixY[x2][y2];

    }
  }
Info << "File matrix_y opened! Reading data from file into coordsFileY" << endl;
}
else
	{
		Info << "File matrix_y not found! Abort" << endl;
		exit(1);
	}


if(ovalXCoords.is_open())
{

	Info << "File ovalXCoords opened! Reading data from file into ovalXCoords" << endl;
for (int xx1 = 0; xx1 < c; xx1++)
	{

    for (int yy1 = 0; yy1 < d; yy1++)
	{

        ovalXCoords >> matrixOvalX[xx1][yy1];

    }
  }
				}

else
	{
		Info << "File ovalXCoords not found! Abort" << endl;
		exit(1);
	}

if(ovalYCoords.is_open())
{

	//Info << "File " << ovalYCoords.c_str() <<" opened! Reading data from file into matrix_y" << endl;
for (int xx2 = 0; xx2 < c; xx2++)
	{

    for (int yy2 = 0; yy2 < d; yy2++)
	{

        ovalYCoords >> matrixOvalY[xx2][yy2];

    }
  }
Info << "File ovalYCoords opened! Reading data from file into ovalYCoords" << endl;
}
else
	{
		Info << "File matrix_y not found! Abort" << endl;
		exit(1);
	}


if(roundXCoords.is_open())
{

	Info << "File ovalXCoords opened! Reading data from file into roundXCoords" << endl;
for (int xxx1 = 0; xxx1 < c; xxx1++)
	{

    for (int yyy1 = 0; yyy1 < d; yyy1++)
	{

        roundXCoords >> matrixRoundX[xxx1][yyy1];

    }
  }
				}

else
	{
		Info << "File roundXCoords not found! Abort" << endl;
		exit(1);
	}

if(roundYCoords.is_open())
{

	//Info << "File " << ovalYCoords.c_str() <<" opened! Reading data from file into matrix_y" << endl;
for (int xxx2 = 0; xxx2 < c; xxx2++)
	{

    for (int yyy2 = 0; yyy2 < d; yyy2++)
	{

        roundYCoords >> matrixRoundY[xxx2][yyy2];

    }
  }
Info << "File roundYCoords opened! Reading data from file into roundYCoords" << endl;
}
else
	{
		Info << "File roundYCoords not found! Abort" << endl;
		exit(1);
	}




while (intPhase < 100)
	{
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        runTime++;
		intPhase++;
        #include "CourantNo.H"
 		#include "readControls.H"
		#include "setDeltaT.H"
Info << "moving Oval " << endl;

#include "moveOval.H"

Info << "moving BM " << endl;
#include "moveBM.H"




	Info << " The tenth BM Y coordinate for this timestep: " << data_readInY[9] << endl;
//Info << " The first Oval X coordinate for this timestep: " << data_readInX2[0] << " and the 3rd x coordinate: " << data_readInX2[4] <<endl;
//Info << " The first Oval Y coordinate for this timestep: " << data_readInY2[0] << " and the 3rd Y coordinate: " << data_readInY2[4] <<endl;

        mesh.update();

 // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;

        if (mesh.changing() && correctPhi)
        {
            #include "correctPhi.H"
        }

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }


		particles.move(g);

		Info<< "Cloud size= "<< particles.size() <<endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
}

// ************************************************************************* //
