#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <ctime>


#include "Neuron.h"


using namespace std;


Neuron::Neuron(int number, string population, double K,
  double externalRateFactor){
  this->number = number;
  this->population = population;
  this->externalRateFactor = externalRateFactor;

  m_0 = 0.5;
  last_spike_time =0;
  last_update_time = 0;
  current_spike_time =0;
  stateSum=0;
  decay_constant = 1;
  threshold_at_last_spike = 0;
  adapation_jump = 0.3;


  // randomly generate delta in (0,1)
  double r = ((double)rand()/(double)RAND_MAX);
  delta = r;

  if (population == "E"){
    threshold = 1;
    original_threshold = threshold;
    externalInput = externalRateFactor*1*m_0*sqrt(K);
    tau = 1;
  }
  else if (population == "I"){
    threshold = 0.7;
    original_threshold = threshold;
    externalInput = externalRateFactor*0.8*m_0*sqrt(K);
    tau = 0.9;
  }

    state = 0;
    previous_state = 0;
    totalExcitatoryInput =0;
    totalInhibitoryInput =0;
    totalInput = 0;
    time_to_be_updated = (0+delta)*tau;
    EI_Ratio = 0;
    update_count =0;
}

void Neuron::addInputData(pair<double,double> totalInput,
  pair<double,double> excitatoryInput, pair<double,double> inhibitoryInput){
  totalInput_timeSeries.push_back(totalInput);
  excitatoryInput_timeSeries.push_back(excitatoryInput);
  inhibitoryInput_timeSeries.push_back(inhibitoryInput);
}

void Neuron::update_time_to_be_updated(){
  time_to_be_updated+=tau;
}

void Neuron::fire(){
  state = 1;
}

void Neuron::rest(){
  state = 0;
}

void Neuron::updateThresholdSmooth(double currentTime){
  if (last_spike_time ==0){
    return;
  }
  //if the neuron happened to fire during this update
  //then raise the threshold by a constant
  else if (currentTime == last_spike_time){
    if (ISI_data.size()!=0){
      threshold = (threshold_at_last_spike - original_threshold)*
      exp(-decay_constant*ISI_data.back()) +
      original_threshold;
    }
      threshold += adapation_jump;
      threshold_at_last_spike = threshold;
      //record
      thresholdVector.push_back(threshold);
  }
  else if (currentTime != last_spike_time){
    threshold = (threshold_at_last_spike - original_threshold)*
    exp(-decay_constant*(currentTime - last_spike_time))
    + original_threshold;
    //record
    thresholdVector.push_back(threshold);
  }
}

void Neuron::updateThresholdDiscrete(double currentTime, double timeElapsed){
  if (last_spike_time ==0){
    //basically do nothing, because the neuron has
    //not spiked a single time yet
    thresholdVector.push_back(threshold);
  }
  else if (currentTime == last_spike_time){
    threshold = original_threshold + (threshold-
    original_threshold)*exp(-decay_constant*timeElapsed);
    threshold += adapation_jump;
    thresholdVector.push_back(threshold);
  }
  else if (currentTime != last_spike_time){
    threshold = (threshold - original_threshold)*
    exp(-decay_constant*timeElapsed)
    + original_threshold;
    //record
    thresholdVector.push_back(threshold);
  }
}
