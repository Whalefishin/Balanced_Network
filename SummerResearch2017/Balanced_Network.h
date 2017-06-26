#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>

#include "GraphLab/adjacencyListGraph.h"
#include "Neuron.h"
#include "GraphLab/stl/stlPriorityQueue.h"

using namespace std;

class Balanced_Network{

public:
  Balanced_Network(double N_E, double N_I, double K, double externalRateFactor);

  ~Balanced_Network();

  //The accessor methods

  //gets threshold for a single neuron
  double getThreshold(Neuron* n);

  //gets externalInput for a single neuron
  double getExternalInput(Neuron* n);

  //gets total input for a single neuron
  double getTotalInput(Neuron* n);

  //gets N_E and N_I
  double getN_E();
  double getN_I();
  double getK();

  //gets the Jmatrix
  double** getJmatrix();

  //gets the current time for the network.
  double getTime();

  vector<double> getEI_Ratios();
  vector<pair<double,double>> getExcMeanAtv();
  vector<pair<double,double>> getInhMeanAtv();
  pair<double,double> getEM_data_exc();
  pair<double,double> getEM_data_inh();

  pair<double,double> getEM_data_exc2();
  pair<double,double> getEM_data_inh2();


  //gets the time vector, the collection of times for which
  //a neuron is updated.
  //Not actually used in the main.
  //vector<double> getTime_Vector();

  //gets the PQ for the update method.
  //For debugging purpose mostly, can saftely ignore.
  STLPriorityQueue<double,Neuron*>* getPQ();


  //The Neuron methods

  //adds N_E excitatory neurons and N_I inhibitory neurons
  //into the network.
  void addNeurons(int N_E, int N_I);

  //check if a neuron is in the network, returns true iff it is.
  bool containsNeuron(Neuron* n);

  //choose a random neuron in the network.
  Neuron* chooseRandomNeuron();

  //get all the neurons in the network and put them into a
  //vector
  vector<Neuron*> getNeurons();



  //The Connection methods

  //insert a connection between two neurons.
  //paramaters: the source and target neuron, the label (can just be
  //empty), and the connection strength.
  void insertConnection(Neuron* source, Neuron* target, string label, double strength);

  //remove a connection between two neurons
  void removeConnection(Neuron* source, Neuron* target);

  //check if two neurons are connected.
  //returns true iff there is a connection established
  //from the source neuron onto the target neuron.
  bool isConnected(Neuron* source, Neuron* target);

  //get all the incoming connections to a neuron (all the
  //connections with that neuron as the target)
  vector<Edge<Neuron*,string,double> > getIncomingConnections(Neuron* n);


  //Initialization methods

  //randomly generate the J matrix with connectivity index K
  void initializeJmatrix(double J_EE, double J_EI, double J_IE, double J_II);

  //randomly set the states of all the neurons in the network
  //to either 1 or 0
  void initializeNeurons();

  //establish connections between neurons based on the J matrix.
  //on average, each neuron should have incoming connections
  //from K excitatory and K inhibitory neurons.
  void establishConnections();

  //The Heaviside function, returns 1 if the input value >0,
  //0 otherwise
  double Heaviside(double value);

  //Update the network once
  //We can change the neuron_to_record to a vector of neurons
  void update(Neuron* neuron_to_record);

  //The same update function, but only records the data
  //for the EM plot.
  void update2();

  void update3();

  void addEI_Ratios();

  int count =0;

  vector<Neuron*> neuron_Vector;


private:
  //helper methods
  void insertNeuron(Neuron* n);
  void removeNeuron(Neuron* n);
  void record(pair<double,double> totalInput_data, pair<double,double> excitatoryInput_data,
    pair<double,double> inhibitoryInput_data, Neuron* neuron_to_record);

  //deletor methods
  void deleteNeurons();
  void deleteJmatrix();

  Graph<Neuron*,string,double>* network;
  double** Jmatrix;


  //The number of excitatory and inhibitory neurons
  double N_E;
  double N_I;

  double K;
  double externalRateFactor;

  double time;
  //vector<double> time_vector;
  STLPriorityQueue<double,Neuron*>* minimum_time_queue;



  //EI ratio and ISI recording
  vector<double> EI_Ratio_Collection;
  vector<double> network_ISI;
  vector<pair<double,double>> excitatoryActivityTimeSeries;
  vector<pair<double,double>> inhibitoryActivityTimeSeries;

  double excStateSum=0;
  double inhStateSum=0;
  double excUpdateCount =0;
  double inhUpdateCount =0;
  int num_update_done;
};
