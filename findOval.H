	//declare motionPtr
	autoPtr<motionSolver> motionPtr2 = motionSolver::New(mesh);

	//identify oval on mesh
	word patchName_Oval = "ovalWindow";
	label patchID_Oval = mesh.boundaryMesh().findPatchID(patchName_Oval);

	// identify leftmost index of Hensen-cell patch
	labelList labelPatchFound2( mesh.boundaryMesh()[patchID_Oval].meshPoints() );
 	pointField meshPoints2(mesh.points());
	vector coord2 = meshPoints2[labelPatchFound2[0]];

	// initialization for min/max-coordinate search
	double minXCoord2 = coord2[0];
	double minYCoord2 = coord2[1];
	double maxYCoord2 = coord2[1];

	int HCLeftEdgeIndex2;
	int HCRightEdgeIndex2;
	//int HCTopIndex2;
	int HCmaxIndex2;
	int ii2 = 1;
	int jj2 = 1;
	int kk2 = 1;
	forAll(labelPatchFound2, index2)
	{
		vector coord2_ = meshPoints2[labelPatchFound2[index2]];
		if ( index2 % 2 == 0 )
		{
			//Info << index2 << " " << coord2_[0] << " " << coord2_[1] << endl;
		}
		if (coord2_[0] < minXCoord2)
		{
			ii2++;
			minXCoord2 = coord2_[0];
			//Info << ii2 << "minXCoord2 = " << minXCoord2 << endl;
			HCLeftEdgeIndex2 = index2;
		}
		if (coord2_[1] < minYCoord2)
		{
			jj2++;
			minYCoord2 = coord2_[1];
			//Info << jj2 << "minYCoord2 = " << minYCoord2 << endl;
			HCRightEdgeIndex2 = index2;
		}
		if (coord2_[1] > maxYCoord2)
		{
			kk2++;
			maxYCoord2 = coord2_[1];
			//Info << kk2 << "maxYCoord2 = " << maxYCoord2 << endl;
			//HCTopIndex = index2;
		}
		HCmaxIndex2 = index2;
	}



	vector leftEdgeCoord2 = meshPoints2[labelPatchFound2[HCLeftEdgeIndex2]];
	vector rightEdgeCoord2 = meshPoints2[labelPatchFound2[HCRightEdgeIndex2]];
	vector LtoR2 = leftEdgeCoord2 - rightEdgeCoord2;
	Info << "LtoR2 " << LtoR2[0] << " " << LtoR2[1] << " " << LtoR2[2] << endl;
	//int nCells2 = int((HCmaxIndex2+1)/2 + 1);
	double LtoRdist2[HCmaxIndex2+1][2];

	//fill with default value
	for (int ii2 = 0; ii2 < HCmaxIndex2+1; ii2++)
	{
		LtoRdist2[ii2][0] = ii2;
		LtoRdist2[ii2][1] = 1e-06;
	}


	// fill LtoRdist with values
	// LtoRdist[:][0] are indices
	// LtoRdist[:][1] are distances
	forAll(labelPatchFound2, index2)
	{
			vector coord2__ = meshPoints2[labelPatchFound2[index2]] - rightEdgeCoord2;
			LtoRdist2[index2][0] = index2;
			//Info << index << endl;
			LtoRdist2[index2][1] = (coord2__[0]*LtoR2[0] + coord2__[1]*LtoR2[1]);
			LtoRdist2[index2][1] *= LtoRdist2[index2][1]; // make sure it's positive
	}

	//now sort LtoRdist (terrible bubble sort, but good enough here ...)
	int tempIndex2;
	double tempValue2;
	bool swap2 = true;
	for (int ii2 = 0; (ii2 < HCmaxIndex2+1) && swap2; ii2++)
	{
		swap2 = false;
		for (int jj2 = 0; jj2 < HCmaxIndex2; jj2++)
		{
		if ( LtoRdist2[jj2][1] > LtoRdist2[jj2+1][1] )
		{
			// swap indices
			tempIndex2 = LtoRdist2[jj2+1][0];
			LtoRdist2[jj2+1][0] = LtoRdist2[jj2][0];
			LtoRdist2[jj2][0] = tempIndex2;

			// swap values
			tempValue2 = LtoRdist2[jj2+1][1];
			LtoRdist2[jj2+1][1] = LtoRdist2[jj2][1];
			LtoRdist2[jj2][1] = tempValue2;

			swap2 = true;
		}
		}
	}

	for (int ii2 = 0; ii2 < HCmaxIndex2+1; ii2++)
	{
		//Info << ii << endl;
		//Info << LtoRdist2[ii2][0] << " " << LtoRdist2[ii2][1] << endl;
	}
