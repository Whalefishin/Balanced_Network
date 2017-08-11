#pragma once

#include <string>
#include <utility>


#include "stdlib.h"
#include "Statistics.h"


using namespace std;

class Neuron{
public:
  Neuron(int number, string population, double K,double m_0,
    double externalRateFactor, double phi, double lambda, double N_E);

  void addInputData(pair<double,double> totalInput,
    pair<double,double> excitatoryInput, pair<double,double> inhibitoryInput);
  inline void update_time_to_be_updated(){
      time_to_be_updated+=tau;
  }
  void fire();
  void rest();
  void updateThresholdSmooth(double currentTime);
  void updateThresholdDiscrete(bool value, double timeElapsed);

  double getMeanTotalInput();
  double getMeanExcInput();
  double getMeanInhInput();

  bool operator ==(const Neuron& n){
    if (number == n.number && population == n.population){
      return true;
    }
    else{
      return false;
    }
  }

  int number;
  int vectorNumber;
  string population;
  double threshold;
  double externalInput;
  double totalInput;
  double totalExcitatoryInput;
  double totalInhibitoryInput;
  double m_0;
  double externalRateFactor;

  double update_count_neuronal;
  double update_count_neuronal_inf;

  double state;
  double previous_state;

  double time_to_be_updated;
  double delta;
  double tau;

  vector<pair<double,double>> totalInput_timeSeries;
  vector<pair<double,double>> excitatoryInput_timeSeries;
  vector<pair<double,double>> inhibitoryInput_timeSeries;

  double totalInput_Sum;
  double totalInput_Sum_inf;
  //double excitatoryInput_Sum;
  //double inhibitoryInput_Sum;

  vector<double> ISI_data;
  vector<double> spike_times;

  //EI Ratios
  double EI_Ratio;
  double EI_Ratio_inf;
  double update_count;
  double update_count_inf;

  //E Inputs and I Inputs
  double E_Input;
  double I_Input;
  double E_Input_inf;
  double I_Input_inf;

  //ISI variables
  double last_spike_time;
  double current_spike_time;

  //Mean activity variables
  double stateSum;


  //Spike Frequency adapation stuff
  double decay_constant;
  double adaptation_jump;
  double original_threshold;

  double threshold_at_last_spike;
  double last_update_time;

  vector<double> thresholdVector;

  //Threshold SD
  double thresholdSum;
  double thresholdSum_inf;
};
