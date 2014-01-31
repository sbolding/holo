//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : Particle1D.cpp
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
//

#include "Particle1D.h"
#include <iostream>
#include "Controller.h"

Particle1D::Particle1D(Mesh* mesh, RNG* rng, string method_str,
	std::vector<CurrentFaceTally*>& current_face_tallies,
	std::vector<CurrentElementTally*>& current_element_tallies,
	std::vector<FluxFaceTally*>& flux_face_tallies,
	std::vector<FluxElementTally*>& flux_element_tallies
	)
{
	_rng = rng;
	_mesh = mesh;
	_position_mfp = -999999999999.; //initialize outside teh domain to check that source samplign works correctly
	_n_elements = _mesh->getNumElems();
	_current_element = 0;	   //particle needs to be somewhere to initialize material properties
	updateElementProperties(); //intiialize to the material properties of the 0-th element, will likely change once sampling occurs, but must initialize
	_method = HoMethods::method_map.at(method_str);
	_is_dead = true;

	if (HoController::PARTICLE_BALANCE)
	{
		_n_abs = 0;
		_n_leak = 0;
		_n_scat = 0;
		_n_terminations=0;
	}

	//store the tallies correctly
	_current_face_tallies = current_face_tallies;
	_current_element_tallies = current_element_tallies;
	_flux_face_tallies = flux_face_tallies;
	_flux_element_tallies = flux_element_tallies;
}

//Sample a path length in cm
inline double Particle1D::samplePathLength()
{
	return -1.*log(_rng->rand_num())*_mfp_tot;
}

//Sample a path length in units of number of MFP, useful for streaming through many cells
inline double Particle1D::samplePathLengthMFP()
{
	return -1.*log(_rng->rand_num());
}


//Update the particle position based on a path length that has been sampled, make sure hasnt streamed out of cell
void Particle1D::streamAcrossGeometry()
{
	//sample a pathlength
	double path_length_mfp = samplePathLengthMFP();

	//determine horizontal displacement
	double displacement_mfp = path_length_mfp*_mu;
	double new_position_mfp = displacement_mfp +_position_mfp;

	//determine if the particle has left the cell, or not
	while ((new_position_mfp < 0.) || (new_position_mfp > _element_width_mfp)) //particle has left the cell
	{
		double path_start = _position_mfp;
		double path_end;

		//Determine the number of mean free paths remaining to stream
		if (_mu >= 0.0) //streaming to the right
		{
			displacement_mfp += (_position_mfp - _element_width_mfp);
			path_end = _element_width_mfp; //leaves to the right
		}
		else //streaming to the left
		{
			displacement_mfp += _position_mfp;
			path_end = 0.0; //leaves to the left
		}
		scoreElementTally(path_start, path_end);
		leaveElement(); //move to the next element

		if (_is_dead) //then particle has leaked, do not score anything else
		{
			return;
		}
		new_position_mfp = _position_mfp + displacement_mfp;
	}

	//Now position is within the current cell
	scoreElementTally(_position_mfp, new_position_mfp); 
	_position_mfp = new_position_mfp;
	
}



//Determine if a scatter or an absorption, and then do teh appropriate behavior after that, depending on the mode
void Particle1D::sampleCollision()
{
	if (_is_dead)
	{
		return;
	}
	if (_method == HoMethods::HOLO_ECMC) //then a pure absorber problem, end the history
	{
		terminateHistory();
		if (HoController::PARTICLE_BALANCE)
		{
			_n_abs++;
		}
	}
	else if ((_method == HoMethods::HOLO_STANDARD_MC) || (_method == HoMethods::STANDARD_MC)) //usual MC sample which event, sample new direction if scattering
	{
		//determine if a scattering event
		if (_rng->rand_num() < _scat_ratio)
		{
			if (HoController::PARTICLE_BALANCE)
			{
				_n_scat++;
			}
			double mu_scat = sampleAngleIsotropic();
			//use angle addition to get the new scattered cosine
			_mu = _mu*mu_scat - sqrt(1. - mu_scat*mu_scat)*sqrt(1. - _mu*_mu);
		}
		else //Non scattering event
		{
			//TODO, if there were something besides absorption possible you would put it here
			terminateHistory();
			if (HoController::PARTICLE_BALANCE) //debug stuff
			{
				_n_abs++;
			}
		}
	}
	else	
	{
		std::cerr << "Input an incorrect mode of operation\n";
		exit(1);
	}

}



//Given the source strength on the two end points (left first), sample from it
void Particle1D::sampleLinDiscontSource(std::vector<double> nodal_values)
{
	//Sample the position based on the nodal values, should write a function to get the area of the source from the element somehow

	if (abs(nodal_values[0] - nodal_values[1]) < 1.E-10) //then constant source, sampling is uniform across the cell
	{
		std::cout << "Probably should change how this check works in sampling source for constant across cell" << std::endl;
		_position_mfp = _rng->rand_num()*_element_width_mfp;
	}
	else //need to sample from lin discontinuous source //THIS ROUTINE WORKS
	{
		double left_hat, right_hat; //Normalized nodal values, such that CDF is normalized
		left_hat = 2.0*nodal_values[0] / (nodal_values[1] + nodal_values[0]);
		right_hat = 2.0 - left_hat;
		//use direct inversion of CDF to sample position, based on quadratic formula
		_position_mfp = -left_hat + sqrt(left_hat*left_hat + 2 * _rng->rand_num()*(right_hat - left_hat));
		_position_mfp /= (right_hat - left_hat);
		_position_mfp *= _element_width_mfp; //convert to mfp
	}
}

void Particle1D::sampleSourceParticle()
{
	//Determine if it is volumetric source, or surface source (depending on the mode you are in, may sample scattering source as well)
	//for now just assume a constant source
	//TODO
	//choose the element you are in
	int elem;
	elem = (int)(_rng->rand_num()*_n_elements);
	_current_element = elem;
	cout << "Current element " << _current_element << endl;
	//Get the position of particle's point of origin
	sampleLinDiscontSource(_mesh->getElement(_current_element)->getExtSourceNodalValues());
	initializeHistory();

	//Assume isotropic source distribution
	_mu = sampleAngleIsotropic();

	//Update particle properties for the new cell
	updateElementProperties();
}

//May not need this function if doesn't do more later
inline void Particle1D::initializeHistory()
{
	_weight = 1.0;
	_is_dead = false;
}

inline double Particle1D::sampleAngleIsotropic()
{
	return _rng->rand_num()*2.0 - 1.;
}

void Particle1D::leaveElement()
{
	scoreFaceTally(); //Score the surface tally, based on the current face
	//Contribute to tallies, determine which element you are entering
	if ( (_mu >= 0.0) && (_current_element < (_n_elements-1)) ) //stream to the cell to the right
	{
		_current_element++; //move to the right one cell	
		updateElementProperties();
		_position_mfp = 0.0; //at the left edge of the new cell
		
	}
	else if( (_mu < 0.0) && (_current_element > 0) ) //stream to the cell to the left
	{
		_current_element--;
		updateElementProperties();
		_position_mfp = _element_width_mfp; //particle is at right edge of the new cell
	}
	else //particle has left the problem domain
	{
		//cout << "I have leaked from element " << _current_element << endl;
		if (HoController::PARTICLE_BALANCE)
		{
			_n_leak++;
		}
		terminateHistory();
	}
}


void Particle1D::runHistory()
{
	cout << "starting history..." << endl;
	//start history
	sampleSourceParticle();

	while (true) //stream the particle until it is absorbed or leaks
	{
		streamAcrossGeometry();
		sampleCollision();
		if (_is_dead)
		{
			break;
		}
	}

}

//return random numbers for use by source or whoever needs one
inline double Particle1D::getRandNum()
{
	return _rng->rand_num();
}

//for when particle has entered a new cell need to update the properties to the current cell
void Particle1D::updateElementProperties()
{
	Element* element = _mesh->getElement(_current_element);
	if (element->getMaterialID() == _element_mat_ID)
	{
		return; //no need to update
	}
	else
	{
		MaterialConstant mat;
		mat = element->getMaterial();

		//Set the material data
		_sigma_tot = mat.getSigmaT();
		_mfp_tot = 1. / _sigma_tot;
		_sigma_scat = mat.getSigmaS();
		_scat_ratio = _sigma_scat/_sigma_tot;
		_sigma_abs = mat.getSigmaA();
		_element_width_mfp = _mfp_tot*element->getElementDimensions()[0];
		_element_mat_ID = element->getMaterialID();
	}

}

void Particle1D::scoreElementTally(double path_start_mfp, double path_end_mfp)
{
	//Score Element tally, need to convert path_length and volume
	//to cm, rather than mfp
	//int _elem_id = _current_element;
	double path_length_cm = abs(path_start_mfp - path_end_mfp)*_mfp_tot/_mu;
	double volume_cm = _element_width_mfp*_mfp_tot;//*1.0cm*1.0cm = h_x(cm^3)
	double normalized_position = 0.5*(path_start_mfp+path_end_mfp)/_element_width_mfp;
	
	//increment tallies
	_current_element_tallies[_current_element]->incrementScore(_weight,
		path_length_cm, _mu, volume_cm, normalized_position); 
	_flux_element_tallies[_current_element]->incrementScore(_weight,
		path_length_cm, _mu, volume_cm, normalized_position);
}

inline void Particle1D::scoreFaceTally()
{
	//The face tally is scored before you have left the cell, so everything is
	//based on the cell you are leaving, not the cell you are entering

	//determine face based on direction
	int face_id = 0; //Leaving to the left
	if (_mu >= 0.0) //Leaving to the right
	{
		face_id = 1;
	}

	//increment the correct tallies
	int face_index = _mesh->getFaceIndex(_current_element,face_id);
	_current_face_tallies[face_index]->incrementScore(_weight, _mu, 1.0); //for one-d, per cm sq
	_flux_face_tallies[face_index]->incrementScore(_weight, _mu, 1.0);
}

inline void Particle1D::terminateHistory()
{
	//TODO may need to do other stuff here
	if (HoController::PARTICLE_BALANCE)
	{
		_n_terminations++;
	}
	_is_dead = true;
}

void Particle1D::printParticleBalance(int n_hist)
{
	cout << "---------------------------------------------------\n"
		<< "               Particle balance\n"
		<< "---------------------------------------------------\n\n"
		<< "	 Number Created: " << n_hist << endl
		<< "	Number Absorbed: " << _n_abs << endl
		<< "      Number Leaked: " << _n_leak << endl
		<< "    Number Scatters: " << _n_scat << endl
		<< "  Number Terminated: " << _n_terminations << endl;	
}