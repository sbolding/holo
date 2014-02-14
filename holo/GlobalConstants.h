// This header contains useful Global Constants, such as pi
#ifndef _GLOBALCONSTANT_H
#define _GLOBALCONSTANT_H

#include <string>
#include <map>

namespace GlobalConstants
{
	const long double PI = 3.141592653589793238L;
	//const long double FOUR_PI = 4.0L * PI;
	const long double FOUR_PI = 2.0L;
	const double RELATIVE_TOLERANCE = 1.E-14; //smallest number that can be added to 1 (effectively anything smaller is ignored, assuming everything is scaled to 1)
}

namespace HoConstants
{
	const double COSINE_CUTOFF = 0.01; //below this value, use fixup factor for cosine
	const double COSINE_SUBSTITUTE_VALUE = 0.5*COSINE_CUTOFF; //half of the cutoff value by default
}

namespace HoMethods
{
	const unsigned int HOLO_ECMC = 1;
	const unsigned int HOLO_STANDARD_MC = 2;
	const unsigned int STANDARD_MC = 3;
	const unsigned int STRATIFIED_SAMPLING = 2; //2 for stratified, 1 for regular alias sampling
	const unsigned int STANDARD_SAMPLING = 1;

	//map these methods to their ints
	//either "holo-ecmc", 'holo-standard-mc', or 'standard-mc'
	const std::map<std::string, int> method_map = { { "holo-ecmc", HOLO_ECMC }, { "holo-standard-mc", HOLO_STANDARD_MC },
		{ "standard-mc", STANDARD_MC } };

	const std::map<std::string, unsigned int> sampling_map = { { "stratified", STRATIFIED_SAMPLING }, { "standard", STANDARD_SAMPLING } };
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