// This file contains all the output control setting.  Eventually this should be done
// in a XML format of some kind TODO

#ifndef _CONTROLLER_H
#define _CONTROLLER_H



namespace LoController
{
	const bool WRITE_MATRIX = false;
	const bool WRITE_LOAD_VECTOR = true;
	const bool WRITE_SOLUTION = true;
}

namespace HoController
{
	const bool PARTICLE_BALANCE = false;
	const bool ADAPTIVE_REFINEMENT = false;
	const double FRACTION_CELLS_TO_REFINE = 0.4;

	//output stuff
	const bool WRITE_ALL_ANGULAR_FLUXES = true; //print out angular fluxes after each batch
	const bool WRITE_RESIDUAL_NORMS = true;
	const bool WRITE_HISTORIES_COMPLETE = false;


	//Sampling stuff
	const int INPUT_SEED = 73907;
	const unsigned int SAMPLING_METHOD = 1; //1 is for alias sampling, currently only one implemented
}

#endif