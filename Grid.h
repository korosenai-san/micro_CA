#pragma once
#include <memory>
#include <vector>
#include <random>
#include <map>
#include "Cell.h"
#include <iostream>

enum NUCLEI
{
	RANDOM = 0,
	EVENLY = 1,
	MC = 2
};

enum BOUNDARY
{
	BLANK = 0,
	REFLECTIVE = 1,
	PERIODIC = 2
};

enum NEIGHBOURHOOD
{
	MOORE = 0,
	NEUMANN = 1
};

class Grid
{
public:
	Grid();
	~Grid();
	virtual uint64_t run();
	//virtual void pause();
	virtual void dumpVTK(const char* fname);

protected:
	std::map<uint16_t, uint32_t>  nuclei_count;
	BOUNDARY boundary;
	NEIGHBOURHOOD neighbour;

	//virtual void grid_copy(std::shared_ptr<void>, std::shared_ptr<void>);
	virtual void grid_swap(std::shared_ptr<void>, std::shared_ptr<void>);
	//virtual void populate_nuclei(uint16_t qty, NUCLEI nuclei_distribution);

private:
	
};