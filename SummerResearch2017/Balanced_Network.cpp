#include <string>
#include <iostream>
#include <math.h>
#include <ctime>


#include"Balanced_Network.h"
#include "stdlib.h"

using namespace std;


Balanced_Network::Balanced_Network(double N_E, double N_I, double K,
double externalRateFactor){
//  srand(1);

  time = 0;
  network = new AdjacencyListGraph<Neuron*,string,double>;
  minimum_time_queue = new STLPriorityQueue<double,Neuron*>;

  //adds E excitatory and I inhibitory neurons to the network
  this->N_E = N_E;
  this->N_I = N_I;
  this->K = K;
  this->externalRateFactor = externalRateFactor;

  for (int i=1;i<=N_E;i++){
    Neuron* n = new Neuron(i,"E",K,externalRateFactor);
    network->insertVertex(n);
    neuron_Vector.push_back(n);
  }

  for (int i=1;i<=N_I;i++){
    Neuron* n = new Neuron(i,"I",K,externalRateFactor);
    network->insertVertex(n);
    neuron_Vector.push_back(n);
  }

  //neuron_Vector = network->getVertices();
  //adds all the update times to a priority queue.
  for (int i=0;i<neuron_Vector.size();i++){
    Neuron* toInsert = neuron_Vector[i];
    minimum_time_queue->insert((-1*toInsert->time_to_be_updated),toInsert);
  }

  //initialize the neurons to states 0 or 1
  initializeNeurons();

  double excAtvCount = 0;
  double inhAtvCount = 0;

  for (int i=0;i<neuron_Vector.size();i++){
    Neuron* n = neuron_Vector[i];
    if (n->state ==1){
      if (n->population == "E"){
        excAtvCount = excAtvCount+1;
      }
      else if (n->population == "I"){
        inhAtvCount = inhAtvCount+1;
      }
    }
  }
  //mean activity stuff
  pair<double,double> excData(0,excAtvCount);
  pair<double,double> inhData(0,inhAtvCount);
  excitatoryActivityTimeSeries.push_back(excData);
  inhibitoryActivityTimeSeries.push_back(inhData);
  cout << "Initial active exc: " + to_string(excAtvCount) << endl;
  cout << "Initial active inh: " + to_string(inhAtvCount) << endl;


//neuron_Vector testing stuff
/*
for (int i=0;i<neuron_Vector.size();i++){
  cout << to_string(neuron_Vector[i]->number) +
  neuron_Vector[i]->population << endl;
}
*/

}

Balanced_Network::~Balanced_Network(){
  deleteNeurons();
  deleteJmatrix();
  delete minimum_time_queue;
  delete network;
}

double Balanced_Network::getThreshold(Neuron* n){
  if (network->containsVertex(n)){
    return n->threshold;
  }
  else {
    throw runtime_error("This neuron is not in the network.");
  }
}

double Balanced_Network::getExternalInput(Neuron* n){
  if (network->containsVertex(n)){
    return n->externalInput;
  }
  else {
    throw runtime_error("This neuron is not in the network.");
  }
}

double Balanced_Network::getTotalInput(Neuron* n){
  if (network->containsVertex(n)){
    return n->totalInput;
  }
  else {
    throw runtime_error("This neuron is not in the network.");
  }
}

double Balanced_Network::getN_E(){
  return N_E;
}

double Balanced_Network::getN_I(){
  return N_I;
}

double Balanced_Network::getK(){
  return K;
}

double** Balanced_Network::getJmatrix(){
  return Jmatrix;
}

double Balanced_Network::getTime(){
  return time;
}

vector<double> Balanced_Network::getEI_Ratios(){
  return EI_Ratio_Collection;
}

vector<pair<double,double>> Balanced_Network::getExcMeanAtv(){
  return excitatoryActivityTimeSeries;
}

vector<pair<double,double>> Balanced_Network::getInhMeanAtv(){
  return inhibitoryActivityTimeSeries;
}

pair<double,double> Balanced_Network::getEM_data_exc(){
  //format: external rate vs. Mean atv
  pair<double,double> toReturn;
  //all the neurons in the same network have the same externalInput
  //so we just choose the first one
  toReturn.first = (neuron_Vector[0]->externalRateFactor);
  double neuronSum=0;
  for (int i=0;i<N_E;i++){
    neuronSum += (neuron_Vector[i]->stateSum/neuron_Vector[i]->update_count);
  }
  toReturn.second = neuronSum/(N_E);
  return toReturn;
}

pair<double,double> Balanced_Network::getEM_data_inh(){
  pair<double,double> toReturn;
  toReturn.first = (neuron_Vector[0]->externalRateFactor);
  double neuronSum=0;
  for (int i=N_E;i<neuron_Vector.size();i++){
    neuronSum += (neuron_Vector[i]->stateSum/neuron_Vector[i]->update_count);
  }
  toReturn.second = neuronSum/(N_I);
  return toReturn;
}

pair<double,double> Balanced_Network::getEM_data_exc2(){
  pair<double,double> toReturn;
  toReturn.first = (neuron_Vector[0]->externalRateFactor);
  double size = excitatoryActivityTimeSeries.size();
  for (int i=1600;i< excitatoryActivityTimeSeries.size();i++){
    excStateSum +=excitatoryActivityTimeSeries[i].second;
  }
  cout << excStateSum << endl;
  toReturn.second = excStateSum/N_E/(excUpdateCount-1600);
  return toReturn;
}

pair<double,double> Balanced_Network::getEM_data_inh2(){
  pair<double,double> toReturn;
  toReturn.first = (neuron_Vector[0]->externalRateFactor);
  double size = inhibitoryActivityTimeSeries.size();
  for (int i=400;i< inhibitoryActivityTimeSeries.size();i++){
    inhStateSum += inhibitoryActivityTimeSeries[i].second;
  }
  cout << inhStateSum << endl;
  toReturn.second = inhStateSum/N_I/(inhUpdateCount-400);
  return toReturn;
}

/*
vector<double> Balanced_Network::getTime_Vector(){
  return time_vector;
}
*/

STLPriorityQueue<double,Neuron*>* Balanced_Network::getPQ(){
  return minimum_time_queue;
}

void Balanced_Network::addNeurons(int N_E, int N_I){
  this->N_E = this->N_E + N_E;
  this->N_I = this->N_I + N_I;

  for (int i=1;i<=N_E;i++){
    Neuron* n = new Neuron(i,"E",K,externalRateFactor);
    network->insertVertex(n);
  }

  for (int i=1;i<=N_I;i++){
    Neuron* n = new Neuron(i,"I",K,externalRateFactor);
    network->insertVertex(n);
  }

  neuron_Vector = network->getVertices();
  //adds all the update times to a priority queue.
  for (int i=0;i<neuron_Vector.size();i++){
    Neuron* toInsert = neuron_Vector[i];
    minimum_time_queue->insert((-1*toInsert->time_to_be_updated),toInsert);
  }
}

bool Balanced_Network::containsNeuron(Neuron* n){
  return network->containsVertex(n);
}

Neuron* Balanced_Network::chooseRandomNeuron(){
  if (neuron_Vector.size()==0){
    throw runtime_error("There are no neurons in the network!");
  }
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  int index = r* neuron_Vector.size();
  return neuron_Vector[index];
}

vector<Neuron*> Balanced_Network::getNeurons(){
  return network->getVertices();
}

void Balanced_Network::insertConnection(Neuron* source, Neuron* target, string label, double strength){
  network->insertEdge(source, target, label, strength);
}

void Balanced_Network::removeConnection(Neuron* source, Neuron* target){
  network->removeEdge(source, target);
}

bool Balanced_Network::isConnected(Neuron* source, Neuron* target){
  return network->containsEdge(source, target);
}

vector<Edge<Neuron*,string,double> > Balanced_Network::getIncomingConnections(Neuron* n){
  return network->getIncomingEdges(n);
}


void Balanced_Network::initializeJmatrix(const double J_EE, const double J_EI, const double J_IE, const double J_II){
  //initialize a N_E+N_I square matrix
  //input paramater K is the connectivity index
  Jmatrix = new double*[int (N_E+N_I)];
  for (int i=0;i<N_E+N_I;i++){
    Jmatrix[i] = new double[int(N_E+N_I)];
  }

  for (int i=0;i<N_E+N_I;i++){
    for (int j=0;j<N_E+N_I;j++){
      Jmatrix[i][j] = 0;

      if (i!=j){
      //r is a random float in [0,1]
      float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      //cout << "Jrand is: " + to_string(r) << endl;

      if (j<N_E && i < N_E && r <= K/N_E){
        Jmatrix[i][j] = J_EE/sqrt(K);
      }
      else if (j<N_E && i >= N_E && r <= K/N_E){
        Jmatrix[i][j] = J_IE/sqrt(K);
      }
      else if (j>=N_E && i < N_E && r <= K/N_I){
        Jmatrix[i][j] = J_EI/sqrt(K);
      }
      else if (j>=N_E && i >= N_E && r <= K/N_I){
        Jmatrix[i][j] = J_II/sqrt(K);
      }
    }
    }
  }
}

void Balanced_Network::initializeNeurons(){
  //  srand(1);

  for (int i=0;i<neuron_Vector.size();i++){
    double r=((double)rand()/(double)RAND_MAX);
    //cout << "initializeNeurons r: " + to_string(r) << endl;
    if (r <=.5){
      neuron_Vector[i]->fire();
    }
    else{
      neuron_Vector[i]->rest();
    }
  }
}

void Balanced_Network::establishConnections(){
  if (neuron_Vector.size()==0){
    throw runtime_error("There are currently no neurons in the network.");
  }

  for (int i=0;i<N_E+N_I;i++){
    for (int j=0;j<N_E+N_I;j++){
      if (Jmatrix[i][j]!=0){
        network->insertEdge(neuron_Vector[j],neuron_Vector[i],"default",Jmatrix[i][j]);
      }
    }
  }
}

double Balanced_Network::Heaviside(double value){
  if (value >0){
    return 1;
  }
  else {
    return 0;
  }
}


void Balanced_Network::update(Neuron* neuron_to_record){
  if (neuron_Vector.size()==0){
    throw runtime_error("There are currently no neurons in the network.");
  }

  //for adapatation
  double timeElapsed;

  //choose the minimum one to update(Max of the negatives)
  Neuron* neuron_to_update = minimum_time_queue->removeMax();

  //reset the total input to the neuron being updated to zero.
  neuron_to_update->totalInput = 0;
  neuron_to_update->totalExcitatoryInput =0;
  neuron_to_update->totalInhibitoryInput =0;

  vector<Edge<Neuron*,string,double> > incomingConnections = network->getIncomingEdges(neuron_to_update);
    //the for-loop sums over all the connected neurons
  for (int j=0;j<incomingConnections.size();j++){
    Neuron* incoming_neuron = incomingConnections[j].source;
    double state = incoming_neuron->state;
    double strength = incomingConnections[j].weight;
    //neuron_to_update->totalInput += strength * state;
    if (strength > 0){
      neuron_to_update->totalExcitatoryInput += (strength*state);
    }
    else if (strength <0){
      neuron_to_update->totalInhibitoryInput += (strength*state);
    }
  }

    //below takes care of externalInput and threshold
    neuron_to_update->totalExcitatoryInput += neuron_to_update->externalInput;
    neuron_to_update->totalInput = neuron_to_update->totalExcitatoryInput +
    neuron_to_update->totalInhibitoryInput;

    //the neuron fires if the totalInput is above zero, rests otherwise.
    neuron_to_update->previous_state = neuron_to_update->state;

    neuron_to_update->state =
    Heaviside(neuron_to_update->totalInput- neuron_to_update->threshold);

    //EM plot stuff(sixth plot)
    //neuron_to_update->stateSum += neuron_to_update->state;

    //EI ratio stuff(second plot)
    //the if statement prevents dividing by zero
    if(neuron_to_update->totalInhibitoryInput !=0){
      neuron_to_update->EI_Ratio += (neuron_to_update->totalExcitatoryInput/
        neuron_to_update->totalInhibitoryInput);
    }
    neuron_to_update->update_count++;

    //update the network time
    time = neuron_to_update->time_to_be_updated;
    timeElapsed = time - neuron_to_update->last_update_time;

    //update the neuron time
    neuron_to_update->update_time_to_be_updated();

    //insert the new update time into the priority queue
    //The queue always has the number of elements equal to
    //the number of neurons.
    minimum_time_queue->insert((-1*neuron_to_update->time_to_be_updated),neuron_to_update);


    //ISI stuff
    if (neuron_to_update->previous_state ==0 &&
    neuron_to_update->state ==1){
      if (neuron_to_update->last_spike_time ==0){
        neuron_to_update->last_spike_time = time;
        //record this as a time in which a spike happened
        neuron_to_update->spike_times.push_back(time);
      }
      else {
        neuron_to_update->current_spike_time = time;
        double ISI_time = neuron_to_update->current_spike_time -
        neuron_to_update->last_spike_time;
        neuron_to_update->ISI_data.push_back(ISI_time);
        neuron_to_update->last_spike_time = time;
        //record this as a time in which a spike happened
        neuron_to_update->spike_times.push_back(time);
      }
    }

    //SFA stuff
    neuron_to_update->updateThresholdDiscrete(time, timeElapsed);
    neuron_to_update->last_update_time = time;


    //Mean Activity stuff
    //this says "if the updated neuron went from rest to active"
    if (neuron_to_update->previous_state ==0 &&
    neuron_to_update->state ==1){
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second+1);
        excitatoryActivityTimeSeries.push_back(data);
      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second+1);
        inhibitoryActivityTimeSeries.push_back(data);
      }
    }
    //this says "else if the updated neuron went from active to rest"
    else if (neuron_to_update->previous_state ==1 &&
    neuron_to_update->state ==0){
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second-1);
        excitatoryActivityTimeSeries.push_back(data);
      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second-1);
        inhibitoryActivityTimeSeries.push_back(data);
      }
    }
    else{
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second);
        excitatoryActivityTimeSeries.push_back(data);
      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second);
        inhibitoryActivityTimeSeries.push_back(data);
      }
    }


    if (neuron_to_record == neuron_to_update){
      //record the totalInput for the target neuron only if it's
      //the one being updated
      //data has the format <time, total input>
      pair<double,double> data_total(time,neuron_to_record->totalInput);
      pair<double,double> data_exc(time,neuron_to_record->totalExcitatoryInput);
      pair<double,double> data_inh(time,neuron_to_record->totalInhibitoryInput);

      record(data_total,data_exc,data_inh,neuron_to_record);
    }
}




void Balanced_Network::update2(){
  if (neuron_Vector.size()==0){
    throw runtime_error("There are currently no neurons in the network.");
  }

  //choose the minimum one to update(Max of the negatives)
  Neuron* neuron_to_update = minimum_time_queue->removeMax();

  //reset the total input to the neuron being updated to zero.
  neuron_to_update->totalInput = 0;
  neuron_to_update->totalExcitatoryInput =0;
  neuron_to_update->totalInhibitoryInput =0;

  vector<Edge<Neuron*,string,double> > incomingConnections = network->getIncomingEdges(neuron_to_update);
    //the for-loop sums over all the connected neurons
  for (int j=0;j<incomingConnections.size();j++){
    Neuron* incoming_neuron = incomingConnections[j].source;
    double state = incoming_neuron->state;
    double strength = incomingConnections[j].weight;
    //neuron_to_update->totalInput += strength * state;
    if (strength > 0){
      neuron_to_update->totalExcitatoryInput += (strength*state);
    }
    else if (strength <0){
      neuron_to_update->totalInhibitoryInput += (strength*state);
    }
  }

    //below takes care of externalInput and threshold
    neuron_to_update->totalExcitatoryInput += neuron_to_update->externalInput;
    neuron_to_update->totalInput = neuron_to_update->totalExcitatoryInput +
    neuron_to_update->totalInhibitoryInput;

    //the neuron fires if the totalInput is above zero, rests otherwise.
    neuron_to_update->previous_state = neuron_to_update->state;
    neuron_to_update->state = Heaviside(neuron_to_update->totalInput- neuron_to_update->threshold);

    //EM plot stuff(sixth plot)
    neuron_to_update->stateSum += neuron_to_update->state;
    neuron_to_update->update_count++;

    //update the network time
    time = neuron_to_update->time_to_be_updated;
    //time_vector.push_back(time);

    //update the neuron time
    neuron_to_update->update_time_to_be_updated();

    //insert the new update time into the priority queue
    //The queue always has the number of elements equal to
    //the number of neurons.
    minimum_time_queue->insert((-1*neuron_to_update->time_to_be_updated),neuron_to_update);


}


void Balanced_Network::update3(){
  if (neuron_Vector.size()==0){
    throw runtime_error("There are currently no neurons in the network.");
  }

  //choose the minimum one to update(Max of the negatives)
  Neuron* neuron_to_update = minimum_time_queue->removeMax();

  //reset the total input to the neuron being updated to zero.
  neuron_to_update->totalInput = 0;
  neuron_to_update->totalExcitatoryInput =0;
  neuron_to_update->totalInhibitoryInput =0;

  vector<Edge<Neuron*,string,double> > incomingConnections = network->getIncomingEdges(neuron_to_update);
    //the for-loop sums over all the connected neurons
  for (int j=0;j<incomingConnections.size();j++){
    Neuron* incoming_neuron = incomingConnections[j].source;
    double state = incoming_neuron->state;
    double strength = incomingConnections[j].weight;
    //neuron_to_update->totalInput += strength * state;
    if (strength > 0){
      neuron_to_update->totalExcitatoryInput += (strength*state);
    }
    else if (strength <0){
      neuron_to_update->totalInhibitoryInput += (strength*state);
    }
  }

    //below takes care of externalInput and threshold
    neuron_to_update->totalExcitatoryInput += neuron_to_update->externalInput;
    neuron_to_update->totalInput = neuron_to_update->totalExcitatoryInput +
    neuron_to_update->totalInhibitoryInput;

    //the neuron fires if the totalInput is above zero, rests otherwise.
    neuron_to_update->previous_state = neuron_to_update->state;
    neuron_to_update->state = Heaviside(neuron_to_update->totalInput- neuron_to_update->threshold);


    //update the network time
    time = neuron_to_update->time_to_be_updated;
    //time_vector.push_back(time);

    //update the neuron time
    neuron_to_update->update_time_to_be_updated();

    //insert the new update time into the priority queue
    //The queue always has the number of elements equal to
    //the number of neurons.
    minimum_time_queue->insert((-1*neuron_to_update->time_to_be_updated),neuron_to_update);



    //Mean Activity stuff
    //this says "if the updated neuron went from rest to active"
    if (neuron_to_update->previous_state ==0 &&
    neuron_to_update->state ==1){
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second+1);
        excitatoryActivityTimeSeries.push_back(data);
        excUpdateCount++;


      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second+1);
        inhibitoryActivityTimeSeries.push_back(data);
        inhUpdateCount++;

      }
    }
    //this says "else if the updated neuron went from active to rest"
    else if (neuron_to_update->previous_state ==1 &&
    neuron_to_update->state ==0){
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second-1);
        excitatoryActivityTimeSeries.push_back(data);
        excUpdateCount++;

      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second-1);
        inhibitoryActivityTimeSeries.push_back(data);
        inhUpdateCount++;

      }
    }
    else{
      if (neuron_to_update->population == "E"){
        pair<double,double> data(time,
          excitatoryActivityTimeSeries.back().second);
        excitatoryActivityTimeSeries.push_back(data);
        excUpdateCount++;
      }
      else if (neuron_to_update->population == "I"){
        pair<double,double> data(time,
          inhibitoryActivityTimeSeries.back().second);
        inhibitoryActivityTimeSeries.push_back(data);
        inhUpdateCount++;

      }
    }

}




void Balanced_Network::addEI_Ratios(){
  for (int i=0;i<neuron_Vector.size();i++){
    EI_Ratio_Collection.push_back(neuron_Vector[i]->EI_Ratio /
      neuron_Vector[i]->update_count);
  }
}


//private methods.

void Balanced_Network::insertNeuron(Neuron* n){
  network->insertVertex(n);
}

void Balanced_Network::removeNeuron(Neuron* n){
  network->removeVertex(n);
}

void Balanced_Network::record(pair<double,double> totalInput_data, pair<double,double> excitatoryInput_data,
  pair<double,double> inhibitoryInput_data, Neuron* neuron_to_record){
  neuron_to_record->addInputData(totalInput_data, excitatoryInput_data,inhibitoryInput_data);
}

void Balanced_Network::deleteNeurons(){
  for (int i=0;i<neuron_Vector.size();i++){
    delete neuron_Vector[i];
  }
}

void Balanced_Network::deleteJmatrix(){
  for (int i=0;i<N_I+N_E;i++){
    delete[] Jmatrix[i];
  }
  delete[] Jmatrix;
}