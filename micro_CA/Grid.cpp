#include "Grid.h"

Grid::Grid()
{
}

Grid::~Grid()
{
}

uint64_t Grid::run()
{
	return uint64_t();
}

void Grid::dumpVTK(const char* fname)
{
}

void Grid::grid_swap(std::shared_ptr<void> g1, std::shared_ptr<void> g2)
{
	std::shared_ptr<void> pom = std::move(g1);

	g1 = std::move(g2);
	g2 = std::move(pom);
}
