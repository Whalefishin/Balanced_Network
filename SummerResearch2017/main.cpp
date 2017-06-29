#include <string>
#include <iostream>
#include <utility>
#include <ctime>



#include "Balanced_Network.h"
#include "GraphLab/stl/stlHashTable.h"
//#include "Neuron.h"
#include "stdlib.h"
//#include "Statistics.h"


using namespace std;

int main(){
//srand(time(NULL));
srand(6);

//Various paramaters
double Num_Excitatory_Neurons = 800;
double Num_Inhibitory_Neurons = 200;
const double K = 20;

const double J_EE = 1;
const double J_EI = -2;
const double J_IE = 1;
const double J_II = -1.8;

const int update_steps = 300000;
const int update_times_for_EM = 1;

const double externalRateFactor = 1;

const double adaptation_jump = 0.3;
const double decay_constant = 0.01;


//If you wish to change the neuronal constants, you have to
//go to the constructor in the Neuron.h file

double Num_All_Neurons = Num_Excitatory_Neurons + Num_Inhibitory_Neurons;






Balanced_Network* neural_network =
new Balanced_Network(Num_Excitatory_Neurons,Num_Inhibitory_Neurons, K,
  externalRateFactor,adaptation_jump,decay_constant);

neural_network->initializeJmatrix(J_EE, J_EI, J_IE, J_II);

//neural_network->establishConnections();


double** J = neural_network->getJmatrix();
vector<Neuron*> nVector = neural_network->neuron_Vector;
//choose a random Neuron to record
Neuron* neuron_to_record = neural_network->chooseRandomNeuron();

//manually choose 200I to record
//Neuron* neuron_to_record = neural_network->neuron_Vector.back();
//Neuron* neuron_to_record = nVector[94];

cout<< "The neuron we choose to record is:" +
to_string(neuron_to_record->number) + neuron_to_record->population << endl;






//if we want to update it a set amout of steps

for (int i=0;i<update_steps;i++){
  neural_network->update(neuron_to_record);
}


//if we want to update it till a certain time spot
double stop_time = 10;


//Print out the data in the terminal


vector<pair<double,double> > totTimeSeries = neuron_to_record->totalInput_timeSeries;
vector<pair<double,double> > excTimeSeries = neuron_to_record->excitatoryInput_timeSeries;
vector<pair<double,double> > inhTimeSeries = neuron_to_record->inhibitoryInput_timeSeries;
/*
for (int i=0;i<totTimeSeries.size();i++){
  cout << "Time: " + to_string(totTimeSeries[i].first) +
  " Total Input: " + to_string(totTimeSeries[i].second) +
  " Total ExcInput: " + to_string(excTimeSeries[i].second) +
  " Total InhInput: " + to_string(inhTimeSeries[i].second) << endl;
}
*/


//print E_I ratios


//add the EI ratio data.
neural_network->addEI_Ratios();

vector<double> EI_ratios = neural_network->getEI_Ratios();
/*
for (int i=0;i<EI_ratios.size();i++){
  cout << "The " + to_string(i+1) + "th EI_ratio is: " +
  to_string(EI_ratios[i]) << endl;
}
*/


//print ISI
/*
for (int i=0;i<nVector.size();i++){
  cout << "The " + to_string(i+1) + "th CV is: " +
  to_string(CoefficientVariation(nVector[i]->ISI_data)) << endl;
}
*/


//print mean Activity
vector<pair<double,double>> exc_mean = neural_network->getExcMeanAtv();
vector<pair<double,double>> inh_mean = neural_network->getInhMeanAtv();


//Mean threshold

vector<pair<double,double>> exc_mean_threshold = neural_network->getMeanExcThresholdTimeSeries();
vector<pair<double,double>> inh_mean_threshold = neural_network->getMeanInhThresholdTimeSeries();



//Adaptation Gain - Lamba vs. EI ratios
vector<double> lambaVector;
vector<double> meanEIratioVector;


//Adaptation Gain - Phi vs. EI ratios
vector<double> phiVector;


//External Input vs. Mean Atv


vector<pair<double,double>> EM_dataVector_excitatory;
vector<pair<double,double>> EM_dataVector_inhibitory;

vector<Balanced_Network*> Collection;

for (int i=1;i<=10;i++){
  srand(6);
  // double r=((double)rand()/(double)RAND_MAX);
  Balanced_Network* toInsert = new
  Balanced_Network(Num_Excitatory_Neurons,Num_Inhibitory_Neurons, K,
    externalRateFactor, adaptation_jump,i*0.01);
  toInsert->initializeJmatrix(J_EE, J_EI, J_IE, J_II);
  //toInsert->establishConnections();
  Collection.push_back(toInsert);
}

//the first way of computing the mean activity
/*
for (int i=0;i<10;i++){
  for (int j=0;j<update_times_for_EM;j++){
    Collection[i]->update2();
  }
  pair<double,double> dataPointExc = Collection[i]->getEM_data_exc();
  pair<double,double> dataPointInh = Collection[i]->getEM_data_inh();
  EM_dataVector_excitatory.push_back(dataPointExc);
  EM_dataVector_inhibitory.push_back(dataPointInh);
}
*/


//the second way of computing the mean activity
for (int i=0;i<10;i++){
  for (int j=0;j<update_times_for_EM;j++){
    //Population gain
    //Collection[i]->update3();

    //Adaptation Gain - Lamba vs. EI ratios
    Neuron* temp = Collection[i]->chooseRandomNeuron();
    Collection[i]->update(temp);
  }
  //pop gain
  /*
  pair<double,double> dataPointExc = Collection[i]->getEM_data_exc2();
  pair<double,double> dataPointInh = Collection[i]->getEM_data_inh2();
  EM_dataVector_excitatory.push_back(dataPointExc);
  EM_dataVector_inhibitory.push_back(dataPointInh);
  */

  //Adaptation Gain - Lamba vs. EI ratios

  Collection[i]->addEI_Ratios();
  lambaVector.push_back(Collection[i]->lamba);
  meanEIratioVector.push_back(mean(Collection[i]->getEI_Ratios()));



  //Adaptation Gain - Phi vs. EI ratios
  /*
  Collection[i]->addEI_Ratios();
  phiVector.push_back(Collection[i]->phi);
  meanEIratioVector.push_back(mean(Collection[i]->getEI_Ratios()));
  */
}


cout << "phiVector size is: " + to_string(phiVector.size()) << endl;

//Outputting.

ofstream timeSeriesTxt1("Data/Time.txt");
ofstream timeSeriesTxt2("Data/TotalInput.txt");
ofstream timeSeriesTxt3("Data/ExcitatoryInput.txt");
ofstream timeSeriesTxt4("Data/InhibitoryInput.txt");

ofstream EIRatiosTxt("Data/EI_ratios.txt");

ofstream CVTxt("Data/ISI_CV.txt");

ofstream excMean_ActivityTxt("Data/Exc_Mean_Activity.txt");
ofstream excTimeTxt("Data/Exc_Time.txt");
ofstream inhMean_ActivityTxt("Data/Inh_Mean_Activity.txt");
ofstream inhTimeTxt("Data/Inh_Time.txt");

ofstream excEMAtvTxt("Data/Exc_EM_activity.txt");
ofstream excEMRateTxt("Data/Exc_EM_external_rate.txt");
ofstream inhEMAtvtxt("Data/Inh_EM_activity.txt");
ofstream inhEMRateTxt("Data/Inh_EM_external_rate.txt");

ofstream thresholdTxt("Data/Threshold.txt");
ofstream spikeTimesTxt("Data/Spike_times.txt");

ofstream parametersTxt("Data/Parameters.txt");

ofstream excMeanThresholdTxt("Data/Exc_Mean_Threshold.txt");
ofstream excThresTimeTxt("Data/Exc_Thres_Time.txt");
ofstream inhMeanThresholdTxt("Data/Inh_Mean_Threshold.txt");
ofstream inhThresTimeTxt("Data/Inh_Thres_Time.txt");

ofstream adaptationGainTxt1("Data/Adaptation_lamba.txt");
ofstream adaptationGainTxt2("Data/Adaptation_EI_Ratios.txt");
ofstream adaptationGainTxt3("Data/Adaptation_Phi.txt");

ofstream plateauTimeTxt("Data/Plateau_Times.txt");
ofstream plateauNeuronsTxt("Data/Plateau_Neurons.txt");


//Inputs
//format: time, total input, excitatory input, inhibitory input
for (int i=0;i<totTimeSeries.size();i++){
  timeSeriesTxt1 <<  to_string(totTimeSeries[i].first) << endl;
  timeSeriesTxt2 << to_string(totTimeSeries[i].second)
<< endl;
timeSeriesTxt3 << to_string(excTimeSeries[i].second)
<< endl;
timeSeriesTxt4 << to_string(inhTimeSeries[i].second) << endl;
}



//EI ratios
for (int i=0;i<EI_ratios.size();i++){
  EIRatiosTxt << to_string(EI_ratios[i]) << endl;
}

//ISI CVs
for (int i=0;i<nVector.size();i++){
  CVTxt << to_string(CoefficientVariation(nVector[i]->ISI_data)) << endl;
}

//Mean activities
for (int i=1000;i<exc_mean.size();i++){
  excTimeTxt << exc_mean[i].first << endl;
  excMean_ActivityTxt << (exc_mean[i].second/Num_Excitatory_Neurons) << endl;
}

for (int i=1000;i<inh_mean.size();i++){
  inhTimeTxt << inh_mean[i].first << endl;
  inhMean_ActivityTxt << (inh_mean[i].second/Num_Inhibitory_Neurons) << endl;
}


//External Input vs. Mean activity
for (int i=0;i<EM_dataVector_excitatory.size();i++){
  excEMRateTxt << EM_dataVector_excitatory[i].first << endl;
  excEMAtvTxt << EM_dataVector_excitatory[i].second << endl;
  inhEMRateTxt << EM_dataVector_inhibitory[i].first << endl;
  inhEMAtvtxt << EM_dataVector_inhibitory[i].second << endl;
}


//Threshold
//each data correspond to the threshold at the times in Time.txt
for (int i=0;i<neuron_to_record->thresholdVector.size();i++){
  thresholdTxt << neuron_to_record->thresholdVector[i] << endl;
}

//spike times
for (int i=0;i<neuron_to_record->spike_times.size();i++){
  spikeTimesTxt << neuron_to_record->spike_times[i] << endl;
}

//Mean Threshold
for (int i=0;i<exc_mean_threshold.size();i++){
  excThresTimeTxt << exc_mean_threshold[i].first << endl;
  excMeanThresholdTxt << exc_mean_threshold[i].second
  /Num_Excitatory_Neurons << endl;
}

for(int i=0;i<inh_mean_threshold.size();i++){
  inhThresTimeTxt << inh_mean_threshold[i].first << endl;
  inhMeanThresholdTxt << inh_mean_threshold[i].second
  /Num_Inhibitory_Neurons << endl;
}


//Adaptation Gain - Lamba vs. EI ratios

for (int i=0;i<lambaVector.size();i++){
  adaptationGainTxt1 << lambaVector[i] << endl;
  adaptationGainTxt2 << meanEIratioVector[i] << endl;
}


//Adaptation Gain - Phi vs. EI ratios
for (int i=0;i<phiVector.size();i++){
  adaptationGainTxt3 << phiVector[i]<< endl;
  adaptationGainTxt2 << meanEIratioVector[i] << endl;
}


//Investigate plateau region
for (int i=0;i<neural_network->activeNeurons.size();i++){
  plateauTimeTxt << neural_network->activeNeurons[i].first << endl;
  plateauNeuronsTxt << to_string(i) << endl;
  plateauNeuronsTxt << " " << endl;
  vector<Neuron*> temp = neural_network->activeNeurons[i].second;
  for (int j=0;j<temp.size();j++){
    plateauNeuronsTxt << to_string(temp[j]->number) + temp[j]->population
    << endl;
  }
}

//parameters
parametersTxt << "This trial run has the following parameters: " << endl;
parametersTxt << "Number of neurons: " + to_string(Num_All_Neurons) << endl;
parametersTxt << "Neuron ratios(Exc/Inh): " + to_string(Num_Excitatory_Neurons/Num_Inhibitory_Neurons) << endl;
parametersTxt << "K: " + to_string(K) << endl;
parametersTxt << "m_0: " + to_string(nVector[0]->m_0) << endl;
parametersTxt << "Number of updates: " + to_string(update_steps) <<endl;
parametersTxt << "Adaptation Jump: " + to_string(nVector[0]->adaptation_jump) << endl;
parametersTxt << "Lamba: " + to_string(nVector[0]->decay_constant) << endl;




//Testing out 95E
/*
vector<Edge<Neuron*,string,double>> incomingConnections =
neural_network->getIncomingConnections(nVector[911]);

cout << "Size is: " + to_string(incomingConnections.size()) << endl;

int inhSum=0;
int excSum=0;

for (int i=0;i<incomingConnections.size();i++){
  Neuron* incoming_neuron = incomingConnections[i].source;
  if (incoming_neuron->population == "I" ){
    inhSum++;
  }
  if (incoming_neuron->population == "E"){
    excSum++;
  }
  cout << to_string(incoming_neuron->number) +
  incoming_neuron->population << endl;
}

cout << "Exc count: " + to_string(excSum) << endl;
cout << "Inh count: " + to_string(inhSum) << endl;

*/


//see if neurons are connected
/*
if (neural_network->isConnected(nVector[94],nVector[760])){
  cout << "95E -> 761E" << endl;
}

if (neural_network->isConnected(nVector[760],nVector[94])){
  cout << "95E <- 761E" << endl;
}

if (neural_network->isConnected(nVector[760],nVector[912])){
  cout << "113I <- 761E" << endl;
}

if (neural_network->isConnected(nVector[912],nVector[760])){
  cout << "113I -> 761E" << endl;
}

if (neural_network->isConnected(nVector[94],nVector[912])){
  cout << "113I <- 95E" << endl;
}

if (neural_network->isConnected(nVector[912],nVector[94])){
  cout << "113I -> 95E" << endl;
}
*/

/*
cout << "One of the entries should be: lamba = "
+ to_string(neural_network->lamba) + "meanEIRatio = " +
to_string(mean(neural_network->getEI_Ratios())) << endl;
*/



//Check if the connections are established appropriately
/*
Neuron* testNeuron = nVector.back();

vector<Edge<Neuron*,string,double>> testEdges =
 neural_network->getIncomingConnections(testNeuron);

 cout << "The test neuron is connected to " +
  to_string(testEdges.size()) + "neurons" << endl;

for(int i=0;i<testEdges.size();i++){
  cout << testEdges[i].weight << endl;
}
*/

//check if the neurons across different networks have the same
//parameters
/*
vector<Neuron*> temp = Collection[0]->neuron_Vector;

for (int i=1;i<10;i++){
  vector<Neuron*> temp2 = Collection[i]->neuron_Vector;
  if (temp.size()!=temp2.size()){
    throw runtime_error("Size not the same");
  }
  for (int j=0;j<temp.size();j++){
    if (temp[j]->number != temp2[j]->number){
      throw runtime_error("name is not equal");
    }
    else if (temp[j]->delta != temp2[j]->delta){
      throw runtime_error("delta is not equal");
    }
    else if (temp[j]->update_count !=temp2[j]->update_count){
      throw runtime_error("update count not equal");
    }
    else {
      cout << "Everything seems fine." << endl;
    }
  }
}
*/

//some ISI testing stuff
/*
for (int j=0;j<nVector.size();j++){
for (int i=0;i<nVector[j]->ISI_data.size();i++){
  cout << "ISI_data for Neuron " + to_string(j+1) + " is: " +
  to_string(nVector[j]->ISI_data[i]) + " ";
}
cout << " " << endl;
}

cout << "Count is: " + to_string(neural_network->count) << endl;
*/

//get all the entires in the Jmatrix
/*
int count = 0;
for (int i=0;i<Num_All_Neurons;i++){
  for (int j=0;j<Num_All_Neurons;j++){
    if (i==Num_All_Neurons-1 && J[i][j]==1){
      count++;
    }
  }
}

cout << count << endl;
*/

//manually check the update times for neurons
/*
for (int i=0;i<nVector.size();i++){
  cout << "The delta for neuron " + to_string(nVector[i]->number)
  + nVector[i]->population + "is:"
  + to_string(nVector[i]->delta) << endl;
}

cout << "The current time is:" + to_string(neural_network->getTime()) << endl;

vector<double> time_vector = neural_network->getTime_Vector();
cout << "Time vector size is:" + to_string(time_vector.size()) << endl;

for (int i=0;i<time_vector.size();i++){
  cout << time_vector[i] << endl;
}
*/

//outputting to a textfile
/*
ofstream A("A.txt");
//write vector entries
for(int i=0;i<Num_All_Neurons;i++)
	{
		for(int j =0; j< Num_All_Neurons;j++)
		{
			A << J[i][j];
			A << " ";
		}
		A << endl;
	}
*/

//test mean and standard deviation and CV
/*
vector<double> data;
data.push_back(10);
data.push_back(5.6);
data.push_back(51);
data.push_back(7.8);
data.push_back(34.5);
data.push_back(17.2);

cout<< "Mean is: " + to_string(mean(data)) << endl;
cout<< "Standard Deviation is: " + to_string(standardDeviation(data)) << endl;
cout << "CV is: " + to_string(CoefficientVariation(data)) << endl;
*/

//check number of excitatory and inhibitory connections
/*
vector<Edge<Neuron*,string,double> > incomingConnections =
neural_network->getIncomingConnections(neuron_to_record);

int excCount =0;
int inhCount = 0;

for (int i=0;i<incomingConnections.size();i++){
  if (incomingConnections[i].weight > 0){
    excCount++;
  }
  else if (incomingConnections[i].weight < 0){
    inhCount++;
  }
}

cout << "Number of excitatory connection is: " +
 to_string(excCount) <<endl;
cout << "Number of inhibitory connection is: " +
 to_string(inhCount) <<endl;
*/



delete neural_network;

for (int i=0;i<Collection.size();i++){
  delete Collection[i];
}

  return 0;
}
