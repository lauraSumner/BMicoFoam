	//declare motionPtr
	autoPtr<motionSolver> motionPtr3 = motionSolver::New(mesh);

	//identify round on mesh
	word patchName_round = "roundWindow";
	label patchID_round = mesh.boundaryMesh().findPatchID(patchName_round);

	// identify leftmost index of round patch
	labelList labelPatchFound3( mesh.boundaryMesh()[patchID_round].meshPoints() );
 	pointField meshPoints3(mesh.points());
	vector coord3 = meshPoints3[labelPatchFound3[0]];

	// initialization for min/max-coordinate search
	double minXCoord3 = coord3[0];
	double minYCoord3 = coord3[1];
	double maxYCoord3 = coord3[1];

	int HCLeftEdgeIndex3;
	int HCRightEdgeIndex3;
	//int HCTopIndex3;
	int HCmaxIndex3;
	int ii3 = 1;
	int jj3 = 1;
	int kk3 = 1;
	forAll(labelPatchFound3, index3)
	{
		vector coord3_ = meshPoints3[labelPatchFound3[index3]];
		if ( index3 % 2 == 0 )
		{
			//Info << index2 << " " << coord2_[0] << " " << coord2_[1] << endl;
		}
		if (coord3_[0] < minXCoord3)
		{
			ii3++;
			minXCoord3 = coord3_[0];
			//Info << ii2 << "minXCoord2 = " << minXCoord2 << endl;
			HCLeftEdgeIndex3 = index3;
		}
		if (coord3_[1] < minYCoord3)
		{
			jj3++;
			minYCoord3 = coord3_[1];
			//Info << jj2 << "minYCoord2 = " << minYCoord2 << endl;
			HCRightEdgeIndex3 = index3;
		}
		if (coord3_[1] > maxYCoord3)
		{
			kk3++;
			maxYCoord3 = coord3_[1];
			//Info << kk2 << "maxYCoord2 = " << maxYCoord2 << endl;
			//HCTopIndex = index2;
		}
		HCmaxIndex3 = index3;
	}



	vector leftEdgeCoord3 = meshPoints3[labelPatchFound3[HCLeftEdgeIndex3]];
	vector rightEdgeCoord3 = meshPoints3[labelPatchFound3[HCRightEdgeIndex3]];
	vector LtoR3 = leftEdgeCoord3 - rightEdgeCoord3;
	Info << "LtoR3 " << LtoR3[0] << " " << LtoR3[1] << " " << LtoR3[2] << endl;
	//int nCells2 = int((HCmaxIndex2+1)/2 + 1);
	double LtoRdist3[HCmaxIndex3+1][2];

	//fill with default value
	for (int ii3 = 0; ii3 < HCmaxIndex3+1; ii3++)
	{
		LtoRdist3[ii3][0] = ii3;
		LtoRdist3[ii3][1] = 1e-06;
	}


	// fill LtoRdist with values
	// LtoRdist[:][0] are indices
	// LtoRdist[:][1] are distances
	forAll(labelPatchFound3, index3)
	{
			vector coord3__ = meshPoints3[labelPatchFound3[index3]] - rightEdgeCoord3;
			LtoRdist3[index3][0] = index3;
			//Info << index << endl;
			LtoRdist3[index3][1] = (coord3__[0]*LtoR3[0] + coord3__[1]*LtoR3[1]);
			LtoRdist3[index3][1] *= LtoRdist3[index3][1]; // make sure it's positive
	}

	//now sort LtoRdist (terrible bubble sort, but good enough here ...)
	int tempIndex3;
	double tempValue3;
	bool swap3 = true;
	for (int ii3 = 0; (ii3 < HCmaxIndex3+1) && swap3; ii3++)
	{
		swap3 = false;
		for (int jj3 = 0; jj3 < HCmaxIndex3; jj3++)
		{
		if ( LtoRdist3[jj3][1] > LtoRdist3[jj3+1][1] )
		{
			// swap indices
			tempIndex3 = LtoRdist3[jj3+1][0];
			LtoRdist3[jj3+1][0] = LtoRdist3[jj3][0];
			LtoRdist3[jj3][0] = tempIndex3;

			// swap values
			tempValue3 = LtoRdist3[jj3+1][1];
			LtoRdist3[jj3+1][1] = LtoRdist3[jj3][1];
			LtoRdist3[jj3][1] = tempValue3;

			swap3 = true;
		}
		}
	}

	for (int ii3= 0; ii3 < HCmaxIndex3+1; ii3++)
	{
		//Info << ii << endl;
		//Info << LtoRdist2[ii2][0] << " " << LtoRdist2[ii2][1] << endl;
	}
