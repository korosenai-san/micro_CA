#include "Grid3D.h"

Grid3D::Grid3D(uint32_t x, uint32_t y, uint32_t z, BOUNDARY bc, NEIGHBOURHOOD nd)
{
	this->x_size = x;
	this->y_size = y;
	this->z_size = z;
	this->boundary = bc;
	this->neighbour = nd;

	grid.reset(new std::shared_ptr<std::shared_ptr<Cell[]>[]>[x_size]);
	grid_Copy.reset(new std::shared_ptr<std::shared_ptr<Cell[]>[]>[x_size]);

	for (size_t i = 0; i < x; i++)
	{
		grid[i].reset(new std::shared_ptr<Cell[]>[y_size]);
		grid_Copy[i].reset(new std::shared_ptr<Cell[]>[y_size]);

		for (size_t j = 0; j < y; j++)
		{
			grid[i][j].reset(new Cell[z_size]());
			grid_Copy[i][j].reset(new Cell[z_size]());
		}
	}
}

Grid3D::~Grid3D()
{
}

uint64_t Grid3D::run()
{
	for (int i = 0; i < x_size; i++)
	{
		for (int j = 0; j < y_size; j++)
		{
#pragma omp parallel for
			for (int k = 0; k < z_size; k++)
			{
				if (grid[i][j][k] == 0)
				{
					std::vector<uint16_t> pom;
					pom = find_neighbours(i, j, k);
					uint16_t neighbour_value;

					switch (neighbour)
					{
					case MOORE:
						neighbour_value = moore(pom);
						if (neighbour_value)
						{
							grid_Copy[i][j][k] = neighbour_value;
#pragma omp critical
							nuclei_count.at(neighbour_value) = nuclei_count.at(neighbour_value) + 1;
						}
						break;
					case NEUMANN:
						if (neumann(pom))
						{
							grid_Copy[i][j][k] = neumann(pom);
#pragma omp critical
							nuclei_count.at(neumann(pom)) = nuclei_count.at(neumann(pom)) + 1;
						}
						break;
					default:
						break;
					}

					pom.clear();
				}
			}
		}
	}

#pragma omp barrier
	grid_copy(grid_Copy, grid);

	uint64_t total_size = 0;

	for (auto entry : nuclei_count)
		total_size += entry.second;

	return x_size * y_size * z_size - total_size;
}

void Grid3D::populate_nuclei(uint16_t qty, NUCLEI nuclei_distribution)
{
	if (nuclei_distribution == RANDOM)
	{
		std::random_device rd;
		std::default_random_engine re(rd());
		std::uniform_int_distribution<uint32_t> uniform_dist_x(0, x_size - 1);
		std::uniform_int_distribution<uint32_t> uniform_dist_y(0, y_size - 1);
		std::uniform_int_distribution<uint32_t> uniform_dist_z(0, z_size - 1);

		for (size_t i = 0; i < qty; i++)
		{
		START_LOOP:

			uint32_t rand_x = uniform_dist_x(re);
			uint32_t rand_y = uniform_dist_y(re);
			uint32_t rand_z = uniform_dist_z(re);

			if (grid[rand_x][rand_y][rand_z] == 0)
			{
				grid[rand_x][rand_y][rand_z] = nuclei_count.size() + 1;
				nuclei_count.insert({ nuclei_count.size() + 1, 1 });
			}
			else
				goto START_LOOP;
		}

		grid_copy(grid, grid_Copy);
	}
		
	if (nuclei_distribution == EVENLY)
	{
		int qty_pom = (int)cbrt(qty);

		for (size_t i = 0; i < qty_pom; i++)
		{
			for (size_t j = 0; j < qty_pom; j++)
			{
				for (size_t k = 0; k < qty_pom; k++)
				{
					grid[(x_size / qty_pom) * i][(y_size / qty_pom) * j][(z_size / qty_pom) * k] = nuclei_count.size() + 1;
					nuclei_count.insert({ nuclei_count.size() + 1, 1 });
				}
			}
		}

		grid_copy(grid, grid_Copy);
	}
		
	if (nuclei_distribution == MC)
	{
		std::random_device rd;
		std::default_random_engine re(rd());
		std::uniform_int_distribution<uint32_t> uniform_dist(1, qty);

		for (size_t i = 0; i < x_size; i++)
		{
			for (size_t j = 0; j < y_size; j++)
			{
				for (size_t k = 0; k < z_size; k++)
				{
					grid[i][j][k] = uniform_dist(re);
					grid_Copy[i][j][k] = 0;
					nuclei_count.insert({ nuclei_count.size() + 1,1 });
				}
			}
		}
	}
}

void Grid3D::run_MC()
{
	std::random_device rd;
	std::default_random_engine re(rd());
	std::uniform_int_distribution<uint32_t> uniform_dist_x(0, x_size - 1);
	std::uniform_int_distribution<uint32_t> uniform_dist_y(0, y_size - 1);
	std::uniform_int_distribution<uint32_t> uniform_dist_z(0, z_size - 1);

	uint64_t iter_limit = x_size * y_size * z_size * 2;
	uint64_t cells_to_visit = x_size * y_size * z_size;

	while (iter_limit && cells_to_visit)
	{
		uint32_t local_x = uniform_dist_x(re);
		uint32_t local_y = uniform_dist_y(re);
		uint32_t local_z = uniform_dist_z(re);

		if (grid_Copy[local_x][local_y][local_z] == 0)
		{
			std::vector<uint16_t> pom;
			pom = find_neighbours(local_x, local_y, local_z);

			int E_1 = [this, pom, local_x, local_y, local_z]() {
				int energy = 0;

				for (auto entry : pom)
					if (this->grid[local_x][local_y][local_z] != entry)
						energy++;

				return energy;
			}();

			std::uniform_int_distribution<uint16_t> uniform_dist_seed(0, pom.size() - 1);
			uint16_t rand_value = uniform_dist_seed(re);

			int E_2 = [pom, rand_value]() {
				int energy = 0;

				for (auto entry : pom)
					if (rand_value != entry)
						energy++;

				return energy;
			}();

			if (E_1 - E_2 <= 0)
				grid[local_x][local_y][local_z] = rand_value;

			grid_Copy[local_x][local_y][local_z] = 1;
			cells_to_visit--;
		}

		iter_limit--;
	}

	for (uint32_t i = 0; i < x_size; i++)
	{
		for (uint32_t j = 0; j < y_size; j++)
		{
			for (uint32_t k = 0; k < z_size; k++)
			{
				if (grid_Copy[i][j][k] == 0)
				{
					std::vector<uint16_t> pom;
					pom = find_neighbours(i, j, k);

					int E_1 = [this, pom, i, j, k]() {
						int energy = 0;

						for (auto entry : pom)
							if (this->grid[i][j][k] != entry)
								energy++;

						return energy;
					}();

					std::uniform_int_distribution<uint16_t> uniform_dist_seed(0, pom.size() - 1);
					uint16_t rand_value = uniform_dist_seed(re);

					int E_2 = [pom, rand_value]() {
						int energy = 0;

						for (auto entry : pom)
							if (rand_value != entry)
								energy++;

						return energy;
					}();

					if (E_1 - E_2 <= 0)
						grid[i][j][k] = rand_value;

					grid_Copy[i][j][k] = 1;
					cells_to_visit--;
				}
			}
		}
	}

	for (uint32_t i = 0; i < x_size; i++)
		for (uint32_t j = 0; j < y_size; j++)
			for (uint32_t k = 0; k < y_size; k++)
				grid_Copy[i][j][k] = 0;

}

void Grid3D::dumpVTK(const char* fname)
{
	FILE* fp;

#pragma warning(suppress : 4996)
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Realy cool data \n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET RECTILINEAR_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", x_size, y_size, z_size);
	fprintf(fp, "X_COORDINATES %d int\n", x_size);
	for (uint32_t i = 0; i < x_size; i++)
		fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fprintf(fp, "Y_COORDINATES %d int\n", y_size);
	for (uint32_t i = 0; i < y_size; i++)
		fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fprintf(fp, "Z_COORDINATES %d int\n", z_size);
	for (uint32_t i = 0; i < z_size; i++)
		fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fprintf(fp, "POINT_DATA %d \n", x_size * y_size * z_size);
	fprintf(fp, "SCALARS data int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (uint32_t i = 0; i < x_size; i++)
		for (uint32_t j = 0; j < y_size; j++)
			for (uint32_t k = 0; k < z_size; k++)
				fprintf(fp, "%d \n", grid[i][j][k].get_value());
	fclose(fp);
}

std::vector<uint16_t> Grid3D::find_neighbours(uint64_t x, uint64_t y, uint64_t z)
{
	std::vector<uint16_t> pom;

	switch (boundary)
	{
	case BLANK:
		if (x == 0 && y == 0 && z == 0)
		{
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][x + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x == 0 && y == 0 && z == z_size - 1)
		{
			pom.push_back(0);	pom.push_back(0);								pom.push_back(0);
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());	pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);	pom.push_back(0);								pom.push_back(0);
			pom.push_back(0);													pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());	pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(0);								pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);								pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);								pom.push_back(0);
		}
		else if (x == x_size - 1 && y == 0 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(0);
		}
		else if (x == x_size - 1 && y == 0 && z == z_size - 1)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][x - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y == y_size - 1 && z == 0)
		{
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(grid[x + 1][y - 1][x + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y == y_size - 1 && z == z_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][x - 1].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
		}
		else if (x == x_size - 1 && y == y_size - 1 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == x_size - 1 && y == y_size - 1 && z == z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y == 0 && z > 0 && z < z_size - 1)
		{
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x == x_size - 1 && y == 0 && z > 0 && z < z_size - 1)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y == 0 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == 0 && z == z_size - 1)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y == y_size - 1 && z > 0 && z < z_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
		}
		else if (x == x_size - 1 && y == y_size - 1 && z > 0 && z < z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1 && z == z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y > 0 && y < y_size - 1 && z == 0)
		{
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(0);
		}
		else if (x == 0 && y > 0 && y < y_size - 1 && z == z_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);	pom.push_back(0);									pom.push_back(0);
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1 && z == z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x == 0 && y > 0 && y < y_size - 1 && z > 0 && z < z < z_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(0);														pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1 && z > 0 && z < z < z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z + 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y == 0 && z > 0 && z < z_size - 1)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1 && z > 0 && z < z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y > 0 && y < y_size - 1 && z == 0)
		{
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y > 0 && y < y_size - 1 && z == z_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
			pom.push_back(0);										pom.push_back(0);									pom.push_back(0);
		}
		else
		{
			pom.push_back(grid[x - 1][y - 1][z - 1].get_value());	pom.push_back(grid[x][y - 1][z - 1].get_value());	pom.push_back(grid[x + 1][y - 1][z - 1].get_value());
			pom.push_back(grid[x - 1][y][z - 1].get_value());		pom.push_back(grid[x][y][z - 1].get_value());		pom.push_back(grid[x + 1][y][z - 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z - 1].get_value());	pom.push_back(grid[x][y + 1][z - 1].get_value());	pom.push_back(grid[x + 1][y + 1][z - 1].get_value());

			pom.push_back(grid[x - 1][y - 1][z].get_value());		pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z].get_value());
			pom.push_back(grid[x - 1][y][z].get_value());																pom.push_back(grid[x + 1][y][z].get_value());
			pom.push_back(grid[x - 1][y + 1][z].get_value());		pom.push_back(grid[x][y + 1][z].get_value());		pom.push_back(grid[x + 1][y + 1][z].get_value());

			pom.push_back(grid[x - 1][y - 1][z + 1].get_value());	pom.push_back(grid[x][y - 1][z].get_value());		pom.push_back(grid[x + 1][y - 1][z + 1].get_value());
			pom.push_back(grid[x - 1][y][z + 1].get_value());		pom.push_back(grid[x][y][z + 1].get_value());		pom.push_back(grid[x + 1][y][z + 1].get_value());
			pom.push_back(grid[x - 1][y + 1][z + 1].get_value());	pom.push_back(grid[x][y + 1][z + 1].get_value());	pom.push_back(grid[x + 1][y + 1][z + 1].get_value());
		}

		break;
	case REFLECTIVE:

		break;
	case PERIODIC:

		break;
	default:
		break;
	}

	return pom;
}

void Grid3D::grid_copy(std::shared_ptr<std::shared_ptr<std::shared_ptr<Cell[]>[]>[]> grid_o, std::shared_ptr<std::shared_ptr<std::shared_ptr<Cell[]>[]>[]> grid_c)
{
	for (size_t i = 0; i < x_size; i++)
	{
		for (size_t j = 0; j < y_size; j++)
		{
			for (size_t k = 0; k < z_size; k++)
			{
				grid_c[i][j][k] = grid_o[i][j][k];
			}
		}
	}
}

uint16_t Grid3D::moore(std::vector<uint16_t> neighbours)
{
	std::map<uint16_t, uint32_t> pom;
	uint16_t max = 0;
	uint32_t max_value = 0;

	for (size_t i = 0; i < 26; i++)
	{
		auto qwe = neighbours.at(i);
		auto search = pom.find(qwe);

		if (search == pom.end())
			pom.insert({ neighbours.at(i), 1 });
		else
			pom.at(search->first) = search->second + 1;
			//pom.insert_or_assign(neighbours.at(i), search->second + 1);
	}

	for (auto entry : pom)
	{
		if (entry.second > max_value && entry.first != 0)
		{
			max = entry.first;
			max_value = entry.second;
		}
	}

	return max;
}

uint16_t Grid3D::neumann(std::vector<uint16_t> neighbours)
{
	std::map<uint16_t, uint32_t> pom;
	uint16_t max = 0;
	uint32_t max_value = 0;

	for (size_t i = 0; i < 26; i++)
	{
		if (i == 4 || i == 10 || i == 12 || i == 13 || i == 15 || i == 21)
		{
			auto qwe = neighbours.at(i);
			auto search = pom.find(qwe);

			if (search == pom.end())
				pom.insert({ neighbours.at(i), 1 });
			else
				pom.insert_or_assign(neighbours.at(i), search->second + 1);
		}
	}

	for (auto entry : pom)
	{
		if (entry.second > max_value && entry.first != 0)
		{
			max = entry.first;
			max_value = entry.second;
		}
	}

	return max;
}
