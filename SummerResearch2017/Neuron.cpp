#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <ctime>


#include "Neuron.h"


using namespace std;


Neuron::Neuron(int number, string population, double K, double m_0,
  double externalRateFactor, double phi, double lambda, double N_E){
  this->number = number;
  this->population = population;
  this->externalRateFactor = externalRateFactor;

  this->m_0 = m_0;
  last_spike_time =0;
  last_update_time = 0;
  current_spike_time =0;
  stateSum=0;
  threshold_at_last_spike = 0;

  //constant adapation_jump
  //adaptation_jump = 0.3;
  //decay_constant = 0.01;

  // randomly generate delta in (0,1)
  double r = ((double)rand()/(double)RAND_MAX);
  delta = r;

  if (population == "E"){
    threshold = 1;
    original_threshold = threshold;
    externalInput = externalRateFactor*1*m_0*sqrt(K);
    //externalInput = externalRateFactor*0.8*m_0*sqrt(K);
    tau = 1;
    vectorNumber = number-1;
    //decay_constant = lambda*tau*2;

  }
  else if (population == "I"){
    threshold = 0.7;
    original_threshold = threshold;
    //externalInput = externalRateFactor*1*m_0*sqrt(K);
    externalInput = externalRateFactor*0.8*m_0*sqrt(K);
    tau = 0.9;
    vectorNumber = number + N_E-1;
    //decay_constant = lambda*tau;
  }
  //jump as a percentage of the original_threshold
  //adaptation_jump = 0.3*threshold;
  //decay_constant = 0.01*tau;


  //get adaptation constants as parameters
  decay_constant = lambda;
  adaptation_jump = phi*threshold;


    state = 0;
    previous_state = 0;
    totalExcitatoryInput =0;
    totalInhibitoryInput =0;
    totalInput = 0;
    time_to_be_updated = (0+delta)*tau;

    EI_Ratio = 0;
    EI_Ratio_inf =0;
    update_count =0;
    update_count_inf = 0;
    update_count_neuronal_inf =0;
    update_count_neuronal=0;
    E_Input =0;
    I_Input =0;
    E_Input_inf=0;
    I_Input_inf=0;
    thresholdSum=0;
    thresholdSum_inf=0;
    totalInput_Sum=0;
    totalInput_Sum_inf=0;

}

void Neuron::addInputData(pair<double,double> totalInput,
  pair<double,double> excitatoryInput, pair<double,double> inhibitoryInput){
  totalInput_timeSeries.push_back(totalInput);
  excitatoryInput_timeSeries.push_back(excitatoryInput);
  inhibitoryInput_timeSeries.push_back(inhibitoryInput);
}

/*
void Neuron::update_time_to_be_updated(){
  time_to_be_updated+=tau;
}
*/

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
      threshold += adaptation_jump;
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

void Neuron::updateThresholdDiscrete(bool value,double timeElapsed){
/*
  if (last_spike_time ==0){
    //basically do nothing, because the neuron has
    //not spiked a single time yet
    thresholdVector.push_back(threshold);
  }*/
  if (value == true){
  if (state == 1){
    threshold += adaptation_jump;
    threshold = original_threshold + (threshold-
    original_threshold)*exp(-decay_constant*timeElapsed);
    thresholdVector.push_back(threshold);
  }
  else if (state ==0){
    threshold = (threshold - original_threshold)*
    exp(-decay_constant*timeElapsed)
    + original_threshold;
    //record
    thresholdVector.push_back(threshold);
  }
}
}

double Neuron::getMeanTotalInput(){
  double sum=0;
  double count=0;
  vector<pair<double,double>> vector = totalInput_timeSeries;
  int size = vector.size();
  for (int i=0;i<size;i++){
    sum += vector[i].second;
    count++;
  }
  return sum/count;
}

double Neuron::getMeanExcInput(){
  double sum=0;
  double count=0;
  vector<pair<double,double>> vector = excitatoryInput_timeSeries;
  int size = vector.size();
  for (int i=0;i<size;i++){
    sum += vector[i].second;
    count++;
  }
  return sum/count;
}

double Neuron::getMeanInhInput(){
  double sum=0;
  double count=0;
  vector<pair<double,double>> vector = inhibitoryInput_timeSeries;
  int size = vector.size();
  for (int i=0;i<size;i++){
    sum += vector[i].second;
    count++;
  }
  return sum/count;
}
