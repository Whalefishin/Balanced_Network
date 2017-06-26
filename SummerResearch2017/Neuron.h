#pragma once

#include <string>
#include <utility>


#include "stdlib.h"

using namespace std;

class Neuron{
public:
  Neuron(int number, string population, double K, double externalRateFactor);

  void addInputData(pair<double,double> totalInput,
    pair<double,double> excitatoryInput, pair<double,double> inhibitoryInput);
  void update_time_to_be_updated();
  void fire();
  void rest();
  void updateThresholdSmooth(double currentTime);
  void updateThresholdDiscrete(double currentTime, double timeElapsed);

  bool operator ==(const Neuron& n){
    if (number == n.number && population == n.population){
      return true;
    }
    else{
      return false;
    }
  }

  int number;
  string population;
  double threshold;
  double externalInput;
  double totalInput;
  double totalExcitatoryInput;
  double totalInhibitoryInput;
  double m_0;
  double externalRateFactor;

  double state;
  double previous_state;

  double time_to_be_updated;
  double delta;
  double tau;

  vector<pair<double,double> > totalInput_timeSeries;
  vector<pair<double,double> > excitatoryInput_timeSeries;
  vector<pair<double,double> > inhibitoryInput_timeSeries;

  vector<double> ISI_data;
  vector<double> spike_times;

  double EI_Ratio;
  double update_count;

  //ISI variables
  double last_spike_time;
  double current_spike_time;

  //Mean activity variables
  double stateSum;


  //Spike Frequency adapation stuff
  double decay_constant;
  double adapation_jump;
  double original_threshold;
  //double decay_time_constant;

  double threshold_at_last_spike;
  double last_update_time;

  vector<double> thresholdVector;
};
