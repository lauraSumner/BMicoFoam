


	pointVectorField& prescribedMotion2 = const_cast<pointVectorField&>
					     				(
											mesh.objectRegistry::
											lookupObject<pointVectorField>
											(
											"pointDisplacement"
											)
				             			);


	vectorField& motionUOvalpatch = refCast<vectorField>(prescribedMotion2.boundaryField()[patchID_Oval]);

	int Psize_Oval = motionUOvalpatch.size();
	vectorField dispVals2(Psize_Oval);


	std::ostringstream convert2; // stringstream used for the conversion
	convert2 << intPhase; //add the value of Number to the characters in the stream



	int HCnumberOfCells2 = (HCmaxIndex2+1)/2;
	int n2 = d/((Psize_Oval/2)-1);
	double  L2 = 0.000625;
	double deltaL2 = L2/n2;




	double x_displ2[HCnumberOfCells2];
	double y_displ2[HCnumberOfCells2];

//need an array big enough for the entire text file
	double data_readInX2[n2*HCnumberOfCells2];
	double data_readInY2[n2*HCnumberOfCells2];
	const point& procPoint2 = mesh.points()[0];


//output some stuff to check the values make sense
	Info << "Total Number of Points on oval window: " << Psize_Oval << endl;
	Info << "Number of Cells on Each Processor: "<< HCnumberOfCells2 << endl;
	Info << "Number of Cores:  "<< n2 << endl;
	Info << "HCMaxIndex2: " << HCmaxIndex2 << endl;


//Row of X coordinates for current timestep read in:

//PROBLEM IS HERE!

		for(int i2 = 0; i2 < d; ++i2)
			  {
				  data_readInX2[i2]= matrixOvalX[intPhase-1][i2];
				Info << "data read in for X: " << data_readInX2[i2] << endl;

			  }

//determine the processor currently being worked on, and then assign the appropriate subset of coordinates depending on which processor (and therefore which part of the geometry) (deltaL is the length of the decomposed section so that procPoint[0] is the first x coordinate of the section)
			for(int m2=0; m2<n2;  ++m2)
					{

						if(procPoint2[0]==m2*deltaL2)
						{
							Pout << " x points on processor " << m2 << endl;
							for(int j2 = 0; j2 < (d/n2); ++j2)
								{
									x_displ2[j2]= data_readInX2[j2+ m2*(HCnumberOfCells2)];
									Info << "x_displ2: " << x_displ2[j2] << endl;

								}


						}

					}




//Repeat for Y coordinates:

		for(int k2 = 0; k2 < d; ++k2)
			  {
				  data_readInY2[k2]= matrixOvalY[intPhase-1][k2];

			  }


			for(int m2=0; m2<n2;  ++m2)
					{

						if(procPoint2[0]==m2*deltaL2)
						{
							Pout << " y points on processor " << m2 << endl;
							for(int l2 = 0; l2 < (d/n2); ++l2)
								{
									y_displ2[l2]= data_readInY2[l2+ m2*(HCnumberOfCells2)];
								}


						}


					}




	// set motion along BM patch
	for(int index2 = 0; index2 < HCmaxIndex2+1; index2++)
	{

		if ( index2 % 2 == 0 )
		{
		dispVals2[LtoRdist2[index2][0]].x() = x_displ2[index2/2];
		dispVals2[LtoRdist2[index2+1][0]].x() = x_displ2[index2/2];
		dispVals2[LtoRdist2[index2][0]].y() = y_displ2[index2/2];
		dispVals2[LtoRdist2[index2+1][0]].y() = y_displ2[index2/2];
		dispVals2[LtoRdist2[index2][0]].z() = 0;
		dispVals2[LtoRdist2[index2+1][0]].z() = 0;

		}
		else
			continue;
	}

	prescribedMotion2.boundaryField()[patchID_Oval] == dispVals2;

   mesh.movePoints(motionPtr2->newPoints());

	Info << dispVals2 << endl;




