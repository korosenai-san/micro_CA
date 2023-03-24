#pragma once
#include "Grid.h"

class Grid3D : Grid
{
public:
	Grid3D(uint32_t x, uint32_t y, uint32_t z, BOUNDARY bc, NEIGHBOURHOOD nd);
	~Grid3D();
	uint64_t run() override;
	void populate_nuclei(uint16_t qty, NUCLEI nuclei_distribution);
	void run_MC();
	//void pause();
	void dumpVTK(const char* fname) override;

private:

	uint32_t x_size, y_size, z_size;
	std::shared_ptr<std::shared_ptr<std::shared_ptr<Cell[]>[]>[]> grid, grid_Copy;

	std::vector<uint16_t> find_neighbours(uint64_t x, uint64_t y, uint64_t z);
	void grid_copy(std::shared_ptr<std::shared_ptr<std::shared_ptr<Cell[]>[]>[]> grid_o, std::shared_ptr<std::shared_ptr<std::shared_ptr<Cell[]>[]>[]> grid_c);

	uint16_t moore(std::vector<uint16_t> neighbours);
	uint16_t neumann(std::vector<uint16_t> neighbours);

};
