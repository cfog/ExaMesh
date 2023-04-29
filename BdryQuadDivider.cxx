//  Copyright 2019 by Carl Ollivier-Gooch.  The University of British
//  Columbia disclaims all copyright interest in the software ExaMesh.//
//
//  This file is part of ExaMesh.
//
//  ExaMesh is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as
//  published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  ExaMesh is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with ExaMesh.  If not, see <https://www.gnu.org/licenses/>.

/*
 * BdryQuadDivider.cxx
 *
 *  Created on: July 8, 2019
 *      Author: cfog
 */

#include "BdryQuadDivider.h"

void BdryQuadDivider::divideInterior() {
}

void BdryQuadDivider::createNewCells() {
	// Okay, sure, these aren't actually cells in the usual sense, but so what?
	for (int jj = 0; jj <= nDivs - 1; jj++) {
		int ii = -1;
		for (ii = 0; ii <= nDivs - 1; ii++) {
			emInt vertsNew1[] = { localVerts[ii][jj][0], localVerts[ii + 1][jj][0],
														localVerts[ii + 1][jj + 1][0],
														localVerts[ii][jj + 1][0] };
			m_pMesh->addBdryQuad(vertsNew1);
			assert(
					m_pMesh->numBdryQuads() <= m_pMesh->maxNBdryQuads());
		} // Done with all quads for this row.
	} // Done with this row (constant j)
}
void BdryQuadDivider::setRefinedVerts(QuadFaceVerts &QFV){
	//std::cout<<"My part id is: "<<QFV.getPartid()<<std::endl; 
	for (int ii = 0; ii <= nDivs ; ii++) {
	 	for (int jj = 0; jj <= nDivs; jj++) {
			int trueI;
			int trueJ;
			emInt vert= localVerts[ii][jj][0]; 
			//std::cout<<"Local verts for quads: "<<vert<<std::endl; 
			QFV.setIntVertInd(ii,jj,vert); 
		}
	}
	//std::cout<<"-----------------------------"<<std::endl; 
}
