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
	//Mesh refinement stuff 
	const bool ADAPTIVE_REFINEMENT = true;
	const double FRACTION_CELLS_TO_REFINE = 0.4;
	const bool REFINE_ACROSS_MU_ZERO = false;
	const bool USE_MAX_JUMP_ERROR = true; //use the max jump error on a cell as the jump error

	//output
	const bool PARTICLE_BALANCE = false;
	const bool WRITE_ALL_ANGULAR_FLUXES = true; //print out angular fluxes after each batch
	const bool WRITE_RESIDUAL_NORMS = true;
	const bool WRITE_HISTORIES_COMPLETE = false;

	//Sampling stuff
	const int INPUT_SEED = 73907;
	const unsigned int SAMPLING_METHOD = 1; //1 is for alias sampling, currently only one implemented
}

#endif