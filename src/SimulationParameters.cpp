#include "SimulationParameters.hpp"


int SimulationParameters::index_cell = 0;

int SimulationParameters::getNextIndex(){
  index_cell = index_cell + 1;

  return index_cell;
}
