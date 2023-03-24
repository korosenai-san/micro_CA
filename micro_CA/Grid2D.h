#pragma once
#include "Grid.h"

class Grid2D : Grid
{
public:
	Grid2D(uint32_t x, uint32_t y, BOUNDARY bc, NEIGHBOURHOOD nd);
	~Grid2D();
	uint64_t run() override;
	void populate_nuclei(uint16_t qty, NUCLEI nuclei_distribution);
	void run_MC();
	//void pause();
	void dumpVTK(const char* fname) override;

private:

	uint32_t x_size, y_size;
	std::shared_ptr<std::shared_ptr<Cell[]>[]> grid, grid_Copy;

	std::vector<uint16_t> find_neighbours(uint64_t x, uint64_t y);
	void grid_copy(std::shared_ptr<std::shared_ptr<Cell[]>[]> grid_o, std::shared_ptr<std::shared_ptr<Cell[]>[]> grid_c);

	uint16_t moore(std::vector<uint16_t> neighbours);
	uint16_t neumann(std::vector<uint16_t> neighbours);

};
