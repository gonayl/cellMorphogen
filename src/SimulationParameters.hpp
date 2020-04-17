#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

class SimulationParameters{
public:

  //Simulation
  static constexpr double TIMESTEP = 0.001;
  static constexpr double TIME_OF_SIMULATION = 48;

/*
  //influence the size of the polarization vector
  static constexpr double IMPACT_POLARISATION_EPI_ON_EPI = -0.08;
  static constexpr double IMPACT_POLARISATION_ENDO_ON_EPI = 0.24;
  static constexpr double IMPACT_POLARISATION_LUMEN_ON_EPI = -0.15;
  static constexpr double VEC_POLARISATION_DECREASE = 0.075;
*/

  //Rules for lumen generation
  static constexpr double THRESHOLD_POLARISATION_EPI = 5;
  static constexpr double AGE_DIV_MIN = 0.05;
  static constexpr double AGE_DIV_LUMEN_MIN = 1;
  static constexpr double TIME_BEETWEN_TWO_LUMEN_GENERATION = 0.5;
  static constexpr double NBR_CELL_BETWEEN_TWO_LUMEN = 2;

  //size of the lumen
  static constexpr double LUMEN_SIZE_FACTOR = 0.022;
  static constexpr double AGE_TO_LUMEN_MATURITY = 30;
  static constexpr double SIZE_MIN_LUMEN = 0.2;

  static constexpr bool SIZE_INCREMENTALE = false;
  static constexpr double SURFACE_IMPACT_ON_LUMEN_DERIVATE = 0.002;
  static constexpr double VECTOR_IMPACT_ON_LUMEN_DERIVATE = 0.0002;


  static constexpr bool GENERATE_VESSEL = true;


  //Function generating a new index at each call
  //!!! this is a static function

  static int index_cell;

  static int getNextIndex();


};


#endif
