// This file contains all the output control setting.  Eventually this should be done
// in a XML format of some kind TODO

#ifndef _CONTROLLER_H
#define _CONTROLLER_H



namespace LoController
{
	const bool WRITE_MATRIX = false;
	const bool WRITE_LOAD_VECTOR = false;
	const bool WRITE_SOLUTION = false;
}

namespace HoController
{
	//Mesh refinement stuff 
	const bool ADAPTIVE_REFINEMENT = true;
	const bool FORCE_REFINEMENT = false;
	const double FRACTION_CELLS_TO_REFINE = 0.3;
	const bool REFINE_ACROSS_MU_ZERO = false;
	const bool USE_MAX_JUMP_ERROR = true; //use the max jump error of all sides as error indicator for refinement

	//IF refinement level is not set, max refinement level is ignored
	const bool MAX_REFINEMENT_LEVEL_SET = false;
	const unsigned int MAX_REFINEMENT_LEVEL = 0;

	//output
	const bool PARTICLE_BALANCE = false;
	const bool WRITE_ALL_ANGULAR_FLUXES = false; //print out angular fluxes after each batch
	const bool WRITE_RESIDUAL_NORMS = false;
	const bool WRITE_RELATIVE_ERROR_NORMS = true; //print out the error after each batch, relative to norm of original psi
	const bool WRITE_HISTORIES_COMPLETE = false;
	const bool WRITE_BATCHES_COMPLETE = true;
	const bool WRITE_MESH_EVERY_REFINEMENT = false;
	const bool WRITE_HALF_RANGE_LO_FLUXES = false;
	const bool WRITE_ALL_LO_DATA = false;

	//Sampling stuff
	const int INPUT_SEED = 73907;
	const unsigned int SAMPLING_METHOD = 1; //1 is for alias sampling, currently only one implemented

	//Particle type stuff
	const bool CONT_WGT_DEPOSITION_PARTICLES = false; //CURRENTLY DOESNT WORK, no aborptions take place, the particle weight is just attenuated exponentially
}

#endif