	//declare motionPtr	
	autoPtr<motionSolver> motionPtr = motionSolver::New(mesh);	

	//identify BM on mesh
	word patchName_BM = "BM1";
	label patchID_BM = mesh.boundaryMesh().findPatchID(patchName_BM);

	// identify leftmost index of Hensen-cell patch
	labelList labelPatchFound( mesh.boundaryMesh()[patchID_BM].meshPoints() );
 	pointField meshPoints(mesh.points());
	vector coord = meshPoints[labelPatchFound[0]];
	
	// initialization for min/max-coordinate search
	double minXCoord = coord[0]; 
	double minYCoord = coord[1];
	double maxYCoord = coord[1];

	int HCLeftEdgeIndex;
	int HCRightEdgeIndex;
	//int HCTopIndex;
	int HCmaxIndex;
	int ii = 1;
	int jj = 1;
	int kk = 1;
	forAll(labelPatchFound, index)
	{
		vector coord_ = meshPoints[labelPatchFound[index]];
		if ( index % 2 == 0 )
		{
			//Info << index << " " << coord_[0] << " " << coord_[1] << endl;
		}		
		if (coord_[0] < minXCoord)
		{
			ii++;
			minXCoord = coord_[0];
			//Info << ii << "minXCoord = " << minXCoord << endl;
			HCLeftEdgeIndex = index;
		}
		if (coord_[1] < minYCoord)
		{
			jj++;
			minYCoord = coord_[1];
			//Info << jj << "minYCoord = " << minYCoord << endl;
			HCRightEdgeIndex = index;
		}
		if (coord_[1] > maxYCoord)
		{
			kk++;
			maxYCoord = coord_[1];
			//Info << kk << "maxYCoord = " << maxYCoord << endl;
			//HCTopIndex = index;
		}
		HCmaxIndex = index;
	}

	

	vector leftEdgeCoord = meshPoints[labelPatchFound[HCLeftEdgeIndex]];
	vector rightEdgeCoord = meshPoints[labelPatchFound[HCRightEdgeIndex]];
	vector LtoR = leftEdgeCoord - rightEdgeCoord;
	//Info << "LtoR " << LtoR[0] << " " << LtoR[1] << " " << LtoR[2] << endl;
	//int nCells = int((HCmaxIndex+1)/2 + 1);
	double LtoRdist[HCmaxIndex+1][2];

	//fill with default value
	for (int ii = 0; ii < HCmaxIndex+1; ii++)
	{
		LtoRdist[ii][0] = ii;
		LtoRdist[ii][1] = 1e-06;
	}


	// fill LtoRdist with values
	// LtoRdist[:][0] are indices
	// LtoRdist[:][1] are distances
	forAll(labelPatchFound, index)
	{
			vector coord__ = meshPoints[labelPatchFound[index]] - rightEdgeCoord;
			LtoRdist[index][0] = index; 
			//Info << index << endl;
			LtoRdist[index][1] = (coord__[0]*LtoR[0] + coord__[1]*LtoR[1]);
			LtoRdist[index][1] *= LtoRdist[index][1]; // make sure it's positive
	}

	//now sort LtoRdist (terrible bubble sort, but good enough here ...)
	int tempIndex;
	double tempValue;
	bool swap = true;
	for (int ii = 0; (ii < HCmaxIndex+1) && swap; ii++)
	{
		swap = false;
		for (int jj = 0; jj < HCmaxIndex; jj++)
		{
		if ( LtoRdist[jj][1] > LtoRdist[jj+1][1] )
		{
			// swap indices
			tempIndex = LtoRdist[jj+1][0];
			LtoRdist[jj+1][0] = LtoRdist[jj][0];
			LtoRdist[jj][0] = tempIndex;

			// swap values
			tempValue = LtoRdist[jj+1][1];
			LtoRdist[jj+1][1] = LtoRdist[jj][1];
			LtoRdist[jj][1] = tempValue;

			swap = true;
		}
		}
	}

	for (int ii = 0; ii < HCmaxIndex+1; ii++)
	{
		//Info << ii << endl;
	//	Info << LtoRdist[ii][0] << " " << LtoRdist[ii][1] << endl;
	}
