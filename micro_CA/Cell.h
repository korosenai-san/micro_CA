#pragma once
#include <memory>

class Cell
{
public:
	Cell();
	~Cell();

	void set_value(uint16_t value);
	uint16_t get_value();

	void operator=(uint16_t v);
	bool operator==(uint16_t v);
	bool operator!=(uint16_t v);

private:
	uint16_t value;

};

