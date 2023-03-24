#include "Cell.h"

Cell::Cell()
{
	this->value = 0;
}

Cell::~Cell()
{
}

void Cell::set_value(uint16_t value)
{
	this->value = value;
}

uint16_t Cell::get_value()
{
	return this->value;
}

void Cell::operator=(uint16_t v)
{
	this->value = v;
}

bool Cell::operator==(uint16_t v)
{
	if (this->value == v)
		return true;
	return false;
}

bool Cell::operator!=(uint16_t v)
{
	if (this->value != v)
		return true;
	return false;
}
