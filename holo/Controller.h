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
	const bool PARTICLE_BALANCE = true;
	const int INPUT_SEED = 73907;
}

#endif