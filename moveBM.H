
	

	pointVectorField& prescribedMotion = const_cast<pointVectorField&>
					     				(
											mesh.objectRegistry::
											lookupObject<pointVectorField>
											(
											"pointDisplacement"
											)
				             			);


	vectorField& motionUBMpatch = refCast<vectorField>(prescribedMotion.boundaryField()[patchID_BM]);
	
	int Psize_BM = motionUBMpatch.size();
	vectorField dispVals(Psize_BM);

	
	std::ostringstream convert; // stringstream used for the conversion
	convert << intPhase; //add the value of Number to the characters in the stream



	int HCnumberOfCells = (HCmaxIndex+1)/2;
	int n = b/((Psize_BM/2)-1); 
	double  L = 0.01;
	double deltaL = L/n;
	


  
	double x_displ[HCnumberOfCells];
	double y_displ[HCnumberOfCells];

//need an array big enough for the entire text file
	double data_readInX[n*HCnumberOfCells];
	double data_readInY[n*HCnumberOfCells];
	const point& procPoint = mesh.points()[0];


//output some stuff to check the values make sense
	Info << "Total Number of Points on BM: " << Psize_BM << endl; 
	Info << "Number of Cells on Each Processor: "<< HCnumberOfCells << endl;
	Info << "Number of Cores:  "<< n << endl;


//Row of X coordinates for current timestep read in:

//PROBLEM IS HERE!

		for(int i = 0; i < b; ++i)
			  { 
				  data_readInX[i]= matrixX[intPhase-1][i];


			  }

//determine the processor currently being worked on, and then assign the appropriate subset of coordinates depending on which processor (and therefore which part of the geometry) (deltaL is the length of the decomposed section so that procPoint[0] is the first x coordinate of the section)
			for(int m=0; m<n;  ++m)	
					{	
						
						if(procPoint[0]==m*deltaL)
						{ 
							//Pout << " x points on processor " << m << endl;
							for(int j = 0; j < (b/n); ++j)
								{
									x_displ[j]= data_readInX[j+ m*(HCnumberOfCells)];
							
								}
		 					
		
						}
		
					}
		



//Repeat for Y coordindates:

		for(int k = 0; k < b; ++k)
			  { 
				  data_readInY[k]= matrixY[intPhase-1][k];

			  }


			for(int m=0; m<n;  ++m)	
					{	
						
						if(procPoint[0]==m*deltaL)
						{ 
							//Pout << " x points on processor " << m << endl;
							for(int l = 0; l < (b/n); ++l)
								{
									y_displ[l]= data_readInY[l+ m*(HCnumberOfCells)];
							
								}
		 					
		
						}
		
	
					}




	// set motion along BM patch
	for(int index = 0; index < HCmaxIndex+1; index++)
	{

		if ( index % 2 == 0 )
		{
		dispVals[LtoRdist[index][0]].x() = x_displ[index/2];
		dispVals[LtoRdist[index+1][0]].x() = x_displ[index/2];
		dispVals[LtoRdist[index][0]].y() = y_displ[index/2];
		dispVals[LtoRdist[index+1][0]].y() = y_displ[index/2];
		dispVals[LtoRdist[index][0]].z() = 0;
		dispVals[LtoRdist[index+1][0]].z() = 0;

		}
		else
			continue;
	}

	prescribedMotion.boundaryField()[patchID_BM] == dispVals;

    mesh.movePoints(motionPtr->newPoints());




