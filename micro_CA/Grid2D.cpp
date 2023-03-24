#include "Grid2D.h"

Grid2D::Grid2D(uint32_t x, uint32_t y, BOUNDARY bc, NEIGHBOURHOOD nd)
{
	this->x_size = x;
	this->y_size = y;
	this->boundary = bc;
	this->neighbour = nd;

	grid.reset(new std::shared_ptr<Cell[]>[x_size]);
	grid_Copy.reset(new std::shared_ptr<Cell[]>[x_size]);

	for (size_t i = 0; i < x; i++)
	{
		grid[i].reset(new Cell[y_size]());
		grid_Copy[i].reset(new Cell[y_size]());
	}
}

Grid2D::~Grid2D()
{
}

uint64_t Grid2D::run()
{
	for (uint32_t i = 0; i < x_size; i++)
	{
		for (uint32_t j = 0; j < y_size; j++)
		{
			if (grid[i][j] == 0)
			{
				std::vector<uint16_t> pom;
				pom = find_neighbours(i, j);

				switch (neighbour)
				{
				case MOORE:
					if (moore(pom))
					{
						grid_Copy[i][j] = moore(pom);
						nuclei_count.at(moore(pom)) = nuclei_count.at(moore(pom)) + 1;
					}
					break;
				case NEUMANN:
					if (neumann(pom))
					{
						grid_Copy[i][j] = neumann(pom);
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

	grid_copy(grid_Copy, grid);

	uint64_t total_size = 0;
	for (auto entry : nuclei_count)
		total_size += entry.second;

	return x_size * y_size - total_size;
}

std::vector<uint16_t> Grid2D::find_neighbours(uint64_t x, uint64_t y)
{
	std::vector<uint16_t> pom;

	switch (boundary)
	{
	case BLANK:
		if (x == 0 && y == 0)
		{
			pom.push_back(0);	pom.push_back(0);							pom.push_back(0);
			pom.push_back(0);												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == 0 && y == y_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(0);												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(0);	pom.push_back(0);							pom.push_back(0);
		}
		else if (x == x_size - 1 && y == 0)
		{
			pom.push_back(0);							pom.push_back(0);							pom.push_back(0);
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(0);
		}
		else if (x == x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(0);
			pom.push_back(0);							pom.push_back(0);							pom.push_back(0);
		}
		else if (x == 0 && y > 0 && y < y_size - 1)
		{
			pom.push_back(0);	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(0);												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(0);	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(0);
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(0);
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(0);
		}
		else if (x > 0 && x < x_size - 1 && y == 0)
		{
			pom.push_back(0);							pom.push_back(0);							pom.push_back(0);
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(0);							pom.push_back(0);							pom.push_back(0);
		}
		else
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}

		break;
	case REFLECTIVE:
		if (x == 0 && y == 0)
		{
			pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x][y].get_value());													pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == 0 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x][y].get_value());													pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x + 1][y].get_value());
		}
		else if (x == x_size - 1 && y == 0)
		{
			pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x][y].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());
		}
		else if (x == x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x][y].get_value());
			pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x][y].get_value());
		}
		else if (x == 0 && y > 0 && y < y_size - 1)
		{
			pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x][y].get_value());													pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == 0)
		{
			pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y].get_value());	pom.push_back(grid[x][y].get_value());		pom.push_back(grid[x + 1][y].get_value());
		}
		else
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}

		break;
	case PERIODIC:
		if (x == 0 && y == 0)
		{
			pom.push_back(grid[x_size - 1][y_size - 1].get_value());	pom.push_back(grid[x][y_size - 1].get_value());	pom.push_back(grid[x + 1][y_size - 1].get_value());
			pom.push_back(grid[x_size - 1][y].get_value());															pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x_size - 1][y + 1].get_value());			pom.push_back(grid[x][y + 1].get_value());		pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == 0 && y == y_size - 1)
		{
			pom.push_back(grid[x_size - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x_size - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x_size - 1][0].get_value());	pom.push_back(grid[x][0].get_value());		pom.push_back(grid[x + 1][0].get_value());
		}
		else if (x == x_size - 1 && y == 0)
		{
			pom.push_back(grid[x - 1][y_size - 1].get_value());	pom.push_back(grid[x][y_size - 1].get_value());	pom.push_back(grid[0][y_size - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());														pom.push_back(grid[0][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());		pom.push_back(grid[x][y + 1].get_value());		pom.push_back(grid[0][y + 1].get_value());
		}
		else if (x == x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[0][x - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[0][x].get_value());
			pom.push_back(grid[x - 1][0].get_value());	pom.push_back(grid[x][0].get_value());		pom.push_back(grid[0][0].get_value());
		}
		else if (x == 0 && y > 0 && y < y_size - 1)
		{
			pom.push_back(grid[x_size - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x_size - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x_size - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x == x_size - 1 && y > 0 && y < y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[0][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[0][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[0][y + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == 0)
		{
			pom.push_back(grid[x - 1][y_size - 1].get_value());	pom.push_back(grid[x][y_size - 1].get_value());	pom.push_back(grid[x + 1][y_size - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());														pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());		pom.push_back(grid[x][y + 1].get_value());		pom.push_back(grid[x + 1][y + 1].get_value());
		}
		else if (x > 0 && x < x_size - 1 && y == y_size - 1)
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][0].get_value());	pom.push_back(grid[x][0].get_value());		pom.push_back(grid[x + 1][0].get_value());
		}
		else
		{
			pom.push_back(grid[x - 1][y - 1].get_value());	pom.push_back(grid[x][y - 1].get_value());	pom.push_back(grid[x + 1][y - 1].get_value());
			pom.push_back(grid[x - 1][y].get_value());												pom.push_back(grid[x + 1][y].get_value());
			pom.push_back(grid[x - 1][y + 1].get_value());	pom.push_back(grid[x][y + 1].get_value());	pom.push_back(grid[x + 1][y + 1].get_value());
		}

		break;
	default:
		break;
	}

	return pom;
}

void Grid2D::populate_nuclei(uint16_t qty, NUCLEI nuclei_distribution)
{
	if (nuclei_distribution == RANDOM)
	{
		std::random_device rd;
		std::default_random_engine re(rd());
		std::uniform_int_distribution<uint32_t> uniform_dist_x(0, x_size - 1);
		std::uniform_int_distribution<uint32_t> uniform_dist_y(0, y_size - 1);

		for (size_t i = 0; i < qty; i++)
		{
		START_LOOP:

			uint32_t rand_x = uniform_dist_x(re);
			uint32_t rand_y = uniform_dist_y(re);
			if (grid[rand_x][rand_y] == 0)
			{
				grid[rand_x][rand_y] = nuclei_count.size() + 1;
				nuclei_count.insert({ nuclei_count.size() + 1, 1 });
			}
			else
				goto START_LOOP;
		}

		grid_copy(grid, grid_Copy);
	}

	if (nuclei_distribution == EVENLY)
	{
		int qty_pom = (int)sqrt(qty);

		for (size_t i = 0; i < qty_pom; i++)
		{
			for (size_t j = 0; j < qty_pom; j++)
			{
				grid[(x_size / qty_pom) * i][(y_size / qty_pom) * j] = nuclei_count.size() + 1;
				nuclei_count.insert({ nuclei_count.size() + 1, 1 });
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
				grid[i][j] = uniform_dist(re);
				grid_Copy[i][j] = 0;
				nuclei_count.insert_or_assign(grid[i][j].get_value(), 1);
			}
		}
	}
}

void Grid2D::run_MC()
{
	std::random_device rd;
	std::default_random_engine re(rd());
	std::uniform_int_distribution<uint32_t> uniform_dist_x(0, x_size - 1);
	std::uniform_int_distribution<uint32_t> uniform_dist_y(0, y_size - 1);

	uint64_t iter_limit = x_size * y_size * 2;
	uint64_t cells_to_visit = x_size * y_size;

	while (iter_limit && cells_to_visit)
	{
		uint32_t local_x = uniform_dist_x(re);
		uint32_t local_y = uniform_dist_y(re);

		if (grid_Copy[local_x][local_y] == 0)
		{

			std::vector<uint16_t> pom;
			pom = find_neighbours(local_x, local_y);

			int E_1 = [this, pom, local_x, local_y]() {
				int energy = 0;

				for (auto entry : pom)
					if (this->grid[local_x][local_y] != entry)
						energy++;

				return energy;
			}();

			std::uniform_int_distribution<uint16_t> uniform_dist_seed(0, pom.size() - 1);
			uint16_t rand_value = pom.at(uniform_dist_seed(re));

			int E_2 = [pom, rand_value]() {
				int energy = 0;

				for (auto entry : pom)
					if (rand_value != entry)
						energy++;

				return energy;
			}();

			if (E_1 - E_2 <= 0)
				grid[local_x][local_y] = rand_value;

			grid_Copy[local_x][local_y] = 1;
			cells_to_visit--;
		}

		iter_limit--;
	}

	for (uint32_t i = 0; i < x_size; i++)
	{
		for (uint32_t j = 0; j < y_size; j++)
		{
			if (grid_Copy[i][j] == 0)
			{
				std::vector<uint16_t> pom;
				pom = find_neighbours(i, j);

				int E_1 = [this, pom, i, j]() {
					int energy = 0;

					for (auto entry : pom)
						if (this->grid[i][j] != entry)
							energy++;

					return energy;
				}();

				std::uniform_int_distribution<uint16_t> uniform_dist_seed(0, pom.size() - 1);
				uint16_t rand_value = pom.at(uniform_dist_seed(re));

				int E_2 = [pom, rand_value]() {
					int energy = 0;

					for (auto entry : pom)
						if (rand_value != entry)
							energy++;

					return energy;
				}();

				if (E_1 - E_2 <= 0)
					grid[i][j] = rand_value;

				grid_Copy[i][j] = 1;
				cells_to_visit--;
			}
		}
	}


	for (uint32_t i = 0; i < x_size; i++)
		for (uint32_t j = 0; j < y_size; j++)
			grid_Copy[i][j] = 0;
}

void Grid2D::dumpVTK(const char* fname)
{
	FILE* fp;

#pragma warning(suppress : 4996)
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Realy cool data \n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET RECTILINEAR_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", x_size, y_size, 1);
	fprintf(fp, "X_COORDINATES %d int\n", x_size);
	for (uint32_t i = 0; i < x_size; i++)
		fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fprintf(fp, "Y_COORDINATES %d int\n", y_size);
	for (uint32_t i = 0; i < y_size; i++)
		fprintf(fp, "%d ", i);
	fprintf(fp, "\n");
	fprintf(fp, "Z_COORDINATES 1 int\n");
	fprintf(fp, "0\n");
	fprintf(fp, "POINT_DATA %d \n", x_size * y_size);
	fprintf(fp, "SCALARS data int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (uint32_t i = 0; i < x_size; i++)
		for (uint32_t j = 0; j < y_size; j++)
			fprintf(fp, "%d \n", grid[i][j].get_value());
	fclose(fp);
}

void Grid2D::grid_copy(std::shared_ptr<std::shared_ptr<Cell[]>[]> grid_o, std::shared_ptr<std::shared_ptr<Cell[]>[]> grid_c)
{
	for (size_t i = 0; i < x_size; i++)
	{
		for (size_t j = 0; j < y_size; j++)
		{
			grid_c[i][j] = grid_o[i][j];
		}
	}
}

uint16_t Grid2D::moore(std::vector<uint16_t> neighbours)
{
	std::map<uint16_t, uint32_t> pom;
	uint16_t max = 0;
	uint32_t max_value = 0;

	for (size_t i = 0; i < 8; i++)
	{
		auto qwe = neighbours.at(i);
		auto search = pom.find(qwe);

		if (search == pom.end())
			pom.insert({ neighbours.at(i), 1 });
		else
			pom.insert_or_assign(neighbours.at(i), search->second + 1);
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

uint16_t Grid2D::neumann(std::vector<uint16_t> neighbours)
{
	std::map<uint16_t, uint32_t> pom;
	uint16_t max = 0;
	uint32_t max_value = 0;

	for (size_t i = 0; i < 8; i++)
	{
		if (i == 1 || i == 3 || i == 4 || i == 6)
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
