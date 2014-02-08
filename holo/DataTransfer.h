//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : DataTransfer.h
//  @ Date : 2/7/2014
//  @ Author : 
//
//

#include "Mesh.h"
#include "HoSolver.h"

#if !defined(_DATATRANSFER_H)
#define _DATATRANSFER_H

class DataTransfer
{
protected:
	Mesh* _mesh;
	HoSolver* _ho_solver;
public:
	DataTransfer(HoSolver* HoSolver, Mesh* mesh);
	void updateLoSystem(); //Calculate the new LoData based on the Hodata.
	virtual void calculateLoData(LoData1D & lo_data, int element_id);	//calculate the LoData parameters based on tallies of Ho solver;
};

#endif  //_DATATRANSFER_H
