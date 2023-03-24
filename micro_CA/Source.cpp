#include <iostream>
#include <iomanip>
#include <chrono>
#include "Grid2D.h"
#include "Grid3D.h"

void main()
{
	int dimension = 3;
	int type;
	int neighbourhood;
	int bc;
	int nuclei = 3;
	int nuclei_cout = 12;
	int iter;
	uint32_t x = 100, y = 100, z = 100;

	std::cout << "Input dimension\n1 - 2D\n2 - 3D\n3 - run prepared CA sim\n4 - run prepared MC sim\n";
	//std::cin >> dimension;
	if (dimension == 1 || dimension == 2)
	{
		std::cout << "Input type\n1 - CA\n2 - MC\n";
		std::cin >> type;
		std::cout << "Input neighbourhood type\n1 - Moore\n2 - von Neumann\n";
		std::cin >> neighbourhood;
		std::cout << "Input bounary condition\n1 - nonperidic\n2 - reflective\n3 - periodic\n";
		std::cin >> bc;
		if (type == 1)
		{
			std::cout << "Input nucleation type\n1 - random\n2 - even\n";
			std::cin >> nuclei;
		}
		std::cout << "Input size x, y, (z)\n";
		std::cin >> x;
		std::cin >> y;
		if (dimension == 2)
			std::cin >> z;
		if (type == 2)
		{
			std::cout << "Input number of iterations\n";
			std::cin >> iter;
		}
		std::cout << "Input number of nuclei\n";
		std::cin >> nuclei_cout;
	}

	NUCLEI n = [nuclei]() {
		NUCLEI pom = RANDOM;

		if (nuclei == 1)
			pom = RANDOM;
		if (nuclei == 2)
			pom = EVENLY;
		if (nuclei == 3)
			pom = MC;

		return pom;
	}();

	BOUNDARY b = [bc]() {
		BOUNDARY pom = BLANK;

		if (bc == 1)
			pom = BLANK;
		if (bc == 2)
			pom = REFLECTIVE;
		if (bc == 3)
			pom = PERIODIC;

		return pom;
	}();

	NEIGHBOURHOOD nd = [neighbourhood]() {
		NEIGHBOURHOOD pom = MOORE;

		if (neighbourhood == 1)
			pom = MOORE;
		if (neighbourhood == 2)
			pom = NEUMANN;

		return pom;
	}();

	if (dimension == 1)
	{
		std::cout << "Running...\n";
		if (type == 1)
		{
			Grid2D grid(x, y, b, nd);
			grid.populate_nuclei(nuclei_cout, n);
			grid.dumpVTK("data0.vtk");
			uint16_t is_working = 1;

			do
			{
				static int runs = 1;
				is_working = grid.run();
				char fn[100];
#pragma warning(suppress : 4996)
				sprintf(fn, "data%d.vtk", runs);
				grid.dumpVTK(fn);
				runs++;

			} while (is_working);
		}
		else if (type == 2)
		{
			Grid2D grid(x, y, b, nd);
			grid.populate_nuclei(nuclei_cout, n);
			grid.dumpVTK("data0.vtk");
			uint16_t is_working = 1;

			do
			{
				static int runs = 1;
				grid.run_MC();
				char fn[100];
#pragma warning(suppress : 4996)
				sprintf(fn, "data%d.vtk", runs);
				grid.dumpVTK(fn);
				runs++;

			} while (is_working <= iter);
		}
	}
	else if (dimension == 2)
	{
		std::cout << "Running...\n";
		if (type == 1)
		{
			Grid3D grid(x, y, z, b, nd);
			grid.populate_nuclei(nuclei_cout, n);
			grid.dumpVTK("data0.vtk");
			uint16_t is_working = 1;

			do
			{
				static int runs = 1;
				is_working = grid.run();
				char fn[100];
#pragma warning(suppress : 4996)
				sprintf(fn, "data%d.vtk", runs);
				grid.dumpVTK(fn);
				runs++;

			} while (is_working);
		}
		else if (type == 2)
		{
			Grid3D grid(x, y, z, b, nd);
			grid.populate_nuclei(nuclei_cout, n);
			grid.dumpVTK("data0.vtk");
			uint16_t is_working = 1;

			do
			{
				static int runs = 1;
				grid.run_MC();
				char fn[100];
#pragma warning(suppress : 4996)
				sprintf(fn, "data%d.vtk", runs);
				grid.dumpVTK(fn);
				runs++;

			} while (is_working <= iter);
		}
	}
	else if (dimension == 3)
	{
		int  a = 100, b = 100;
		auto start_time = std::chrono::high_resolution_clock::now();
		std::cout << "Running...\n";
		Grid2D grid(a, b, BLANK, MOORE);
		grid.populate_nuclei(12, EVENLY);
		grid.dumpVTK("data0.vtk");
		uint16_t is_working = 1;

		do
		{
			is_working = grid.run();
			std::cout << "\r[";
			double pom = a * b - is_working;
			for (size_t i = 0; i < 100; i++)
			{
				if (i < (pom/(a*b)*100))
					std::cout << "|";
				else
					std::cout << " ";
			}
			std::cout << "] " << std::setprecision(4) << (pom / (a * b) * 100) << "%";
		} while (is_working);

		grid.dumpVTK("data1.vtk");
		auto end_time = std::chrono::high_resolution_clock::now();
		auto time = end_time - start_time;
		std::cout << "\nReal time: " << time / std::chrono::milliseconds(1) << " ms" << std::endl;
	}
	else if (dimension == 4)
	{
		std::cout << "Running...\n";
		Grid3D grid(50, 50, 50, BLANK, MOORE);
		grid.populate_nuclei(12, MC);
		grid.dumpVTK("data0.vtk");
		uint16_t is_working = 50;

		do
		{
			grid.run_MC();
			is_working--;
		} while (is_working);

		grid.dumpVTK("data1.vtk");
	}

	/*uint64_t pom = 100;

	if (true)
	{
		uint16_t size = pom;
		do
		{
			grid.run_MC();
			char fn[100];
#pragma warning(suppress : 4996)
			sprintf(fn, "data%d.vtk", iter);
			grid.dumpVTK(fn);
			iter++;

			std::cout << "\r[";
			for (size_t i = 0; i < size; i++)
			{
				if (i <= (100 - pom))
					std::cout << "|";
				else
					std::cout << " ";
			}
			std::cout << "]";

			pom--;
		} while (pom);
	}
	else
	{
		do
		{
			pom = grid.run();
			char fn[100];
#pragma warning(suppress : 4996)
			sprintf(fn, "data%d.vtk", iter);
			grid.dumpVTK(fn);
			iter++;


		} while (pom);
	}*/
	
}


