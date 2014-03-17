// Uses alias sampling to sample from a histogram (does not need to be normalized). 
// Needs to be passed the histogram in the form of vector<double>  
// and a random number generator.  Will return the location of bin in the array.
// 
// WARNING: Currently a copy of the array is passed in, be careful if you switch
// to pass by reference as it will destroy the array outside of the class as well
//  Generated by StarUML(tm) C++ Add-In
//
//
//  @ File Name : AliasSampler.h
//  @ Date : 2/14/2013
//  @ Author : Simon R Bolding
//
//


#ifndef _ALIASSAMPLER_H
#define _ALIASSAMPLER_H

#include <iostream>
#include <vector>

struct AliasData
{
	unsigned int alias_event;
	double non_alias_probability;
};

namespace HistogramUtilities
{
	void normalizeDiscreteDistribution(std::vector<double> & bin_probabilities); //will return a normalized pdf, overwrites bin_probabilities
	void scaleDiscreteDistribution(std::vector<double> & bin_probablities, double scalar); //multiply all members of pdf by constant
}

class AliasSampler
{
private:
	std::vector<AliasData> _alias_data; //index of array is the corresponding bin
	void createAliasTable(std::vector<double> &pp);
	int _n_bins;

public:
	AliasSampler(std::vector<double> bin_probabilities, bool normalized=false); 
	AliasSampler(); //Default constructor
	//Need the discrete probabilites of each bin, but can be unnormalized histogram
	
	//Functions
	unsigned int sampleBin(double random_number_0_1, double random_number_2);  //return a sample of which bin to use, based on 2 random numbers between 0 and 1
};

inline unsigned int AliasSampler::sampleBin(double random_number1, double random_number2)
{
	int bin = (int)(random_number1*_n_bins); //which bin are you in
	return (random_number2 < _alias_data[bin].non_alias_probability ? bin : _alias_data[bin].alias_event); //Did you stay in that bin, or go to Alias bin?
}

#endif //_ALIASSAMPLER_H