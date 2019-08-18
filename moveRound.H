


	pointVectorField& prescribedMotion3 = const_cast<pointVectorField&>
					     				(
											mesh.objectRegistry::
											lookupObject<pointVectorField>
											(
											"pointDisplacement"
											)
				             			);


	vectorField& motionUroundpatch = refCast<vectorField>(prescribedMotion3.boundaryField()[patchID_round]);

	int Psize_round = motionUroundpatch.size();
	vectorField dispVals3(Psize_round);


	std::ostringstream convert3; // stringstream used for the conversion
	convert3 << intPhase; //add the value of Number to the characters in the stream



	int HCnumberOfCells3 = (HCmaxIndex3+1)/2;
	int n3 = d/((Psize_round/2)-1);
	double  L3 = 0.000625;
	double deltaL3 = L3/n3;




	double x_displ3[HCnumberOfCells3];
	double y_displ3[HCnumberOfCells3];

//need an array big enough for the entire text file
	double data_readInX3[n3*HCnumberOfCells3];
	double data_readInY3[n3*HCnumberOfCells3];
	const point& procPoint3 = mesh.points()[0];


//output some stuff to check the values make sense
	Info << "Total Number of Points on round window: " << Psize_round << endl;
	Info << "Number of Cells on Each Processor: "<< HCnumberOfCells3 << endl;
	Info << "Number of Cores:  "<< n3 << endl;
	Info << "HCMaxIndex3: " << HCmaxIndex3 << endl;


//Row of X coordinates for current timestep read in:

//PROBLEM IS HERE!

		for(int i3 = 0; i3 < d; ++i3)
			  {
				  data_readInX3[i3]= matrixRoundX[intPhase-1][i3];
				Info << "data read in for X: " << data_readInX3[i3] << endl;

			  }

//determine the processor currently being worked on, and then assign the appropriate subset of coordinates depending on which processor (and therefore which part of the geometry) (deltaL is the length of the decomposed section so that procPoint[0] is the first x coordinate of the section)
			for(int m3=0; m3<n3;  ++m3)
					{

						if(procPoint3[0]==m3*deltaL3)
						{
							Pout << " x points on processor " << m3 << endl;
							for(int j3 = 0; j3 < (d/n3); ++j3)
								{
									x_displ3[j3]= data_readInX3[j3+ m3*(HCnumberOfCells3)];
									Info << "x_displ3: " << x_displ3[j3] << endl;

								}


						}

					}




//Repeat for Y coordinates:

		for(int k3 = 0; k3 < d; ++k3)
			  {
				  data_readInY3[k3]= matrixRoundY[intPhase-1][k3];

			  }


			for(int m3=0; m3<n3;  ++m3)
					{

						if(procPoint3[0]==m3*deltaL3)
						{
							Pout << " y points on processor " << m3 << endl;
							for(int l3 = 0; l3 < (d/n3); ++l3)
								{
									y_displ3[l3]= data_readInY3[l3+ m3*(HCnumberOfCells3)];
								}


						}


					}




	// set motion along BM patch
	for(int index3 = 0; index3 < HCmaxIndex3+1; index3++)
	{

		if ( index3 % 2 == 0 )
		{
		dispVals3[LtoRdist3[index3][0]].x() = x_displ3[index3/2];
		dispVals3[LtoRdist3[index3+1][0]].x() = x_displ3[index3/2];
		dispVals3[LtoRdist3[index3][0]].y() = y_displ3[index3/2];
		dispVals3[LtoRdist3[index3+1][0]].y() = y_displ3[index3/2];
		dispVals3[LtoRdist3[index3][0]].z() = 0;
		dispVals3[LtoRdist3[index3+1][0]].z() = 0;

		}
		else
			continue;
	}

	prescribedMotion3.boundaryField()[patchID_round] == dispVals3;

   mesh.movePoints(motionPtr3->newPoints());

	//Info << dispVals2 << endl;




