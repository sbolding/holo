// This header contains useful Global Constants, such as pi
#ifndef _GLOBALCONSTANT_H
#define _GLOBALCONSTANT_H

#include <string>

namespace GlobalConstants
{
	const long double PI = 3.141592653589793238L;
	//const long double FOUR_PI = 4.0L * PI;
	const long double FOUR_PI = 2.0L;
}

namespace HoMethods
{
	const unsigned int HOLO_ECMC = 0;
	const unsigned int HOLO_STANDARD_MC = 1;
	const unsigned int STANDARD_MC = 2;
}

namespace EquationMaps1D
{
	const unsigned int DOF_MAP[4] = {3 , 2, 0, 1 }; //Map for DOF aliases, this is used in Element1D constructor
	//Order of values is phi_left_plus, phi_right_plus, phi_left_minus, phi_right_minus

	//Map the alias function to actual function, this order is done in the initialization list
	//----- WARNING: DO NOT CHANGE THIS BLOCK -----
	//_phi_left_plus = _elem_dofs[DOF_MAP[0]];
	//_phi_right_plus = _elem_dofs[DOF_MAP[1]];
	//_phi_left_minus = _elem_dofs[DOF_MAP[2]];
	//_phi_right_minus = _elem_dofs[DOF_MAP[3]];
	// --------------------------------------------

	const unsigned int EQN_MAP[4] = { 0, 1, 2, 3 }; //Map for ordering of equations (the row in matrix they go in)
	//Order of values is <.>_left_plus, <.>_right_plus, <.>_left_minus, <.>_right_minus, where each one represents
	//how that moment equation was derived.

	//EQN_MAP[0] refers to < . >_L^+ term
	//EQN_MAP[1] refers to < . >_R^+ term
	//EQN_MAP[2] refers to < . >_L^- term
	//EQN_MAP[3] refers to < . >_R^- term
}

#endif