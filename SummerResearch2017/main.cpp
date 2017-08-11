#include <string>
#include <iostream>
#include <utility>
#include <ctime>



#include "Balanced_Network.h"
#include "GraphLab/stl/stlHashTable.h"
#include "stdlib.h"


using namespace std;

int main(){
//srand(time(NULL));
srand(6);

//Various paramaters
const double Num_Excitatory_Neurons = 8;
const double Num_Inhibitory_Neurons = 2;
const double Num_Excitatory_Neurons_Scale = 8000;
const double Num_Inhibitory_Neurons_Scale = 2000;

const double K = 200;
const double K_Scale = 200;

const double m_0 = 0.2;
const double m_0_Scale = 0.2;

const double J_EE = 1;
const double J_EI = -2;
const double J_IE = 1;
const double J_II = -1.8;

const int update_times = 20;
const int update_times_Scale = 10000*300;

const double Num_Scale1 = 1;
const double Num_Scale2 = 1;

const double externalRateFactor = 1;
const double adaptation_jump = 0;
const double decay_constant = 0.025;

double Num_All_Neurons = Num_Excitatory_Neurons + Num_Inhibitory_Neurons;

//If you wish to change the neuronal constants, you have to
//go to the constructor in the Neuron.h file






Balanced_Network* neural_network =
new Balanced_Network(Num_Excitatory_Neurons,Num_Inhibitory_Neurons, K,
  J_EE, J_EI, J_IE, J_II,m_0,externalRateFactor,adaptation_jump,decay_constant,update_times);

double** J = neural_network->getJmatrix();
vector<Neuron*> nVector = neural_network->neuron_Vector;
//choose a random Neuron to record
Neuron* neuron_to_record = neural_network->chooseRandomNeuron();

cout<< "The neuron we choose to record is:" +
to_string(neuron_to_record->number) + neuron_to_record->population << endl;






//if we want to update it a set amout of steps

for (int i=0;i<update_times;i++){
  //this update collects total input for a single, randomly chosen neuron
  //neural_network->update(neuron_to_record);

  //this update collects total input for each population
  neural_network->update(neuron_to_record);
}


//if we want to update it till a certain time spot
double stop_time = 10;

vector<pair<double,double> > totTimeSeries = neuron_to_record->totalInput_timeSeries;
vector<pair<double,double> > excTimeSeries = neuron_to_record->excitatoryInput_timeSeries;
vector<pair<double,double> > inhTimeSeries = neuron_to_record->inhibitoryInput_timeSeries;

vector<pair<double,double>> totalInputExc = neural_network->totalInput_exc_timeSeries;
vector<pair<double,double>> totalInputInh = neural_network->totalInput_inh_timeSeries;

//add the EI ratio data.
neural_network->addEI_Ratios();

vector<double> EI_ratios = neural_network->getEI_Ratios();

//print mean Activity
vector<pair<double,double>> exc_mean = neural_network->getExcMeanAtv();
vector<pair<double,double>> inh_mean = neural_network->getInhMeanAtv();

//Mean threshold

vector<pair<double,double>> exc_mean_threshold = neural_network->getMeanExcThresholdTimeSeries();
vector<pair<double,double>> inh_mean_threshold = neural_network->getMeanInhThresholdTimeSeries();



//Parameter Scaling plots


//Phi & Lambda -> x & y
vector<double> lambdaVector;
vector<double> phiVector;

//EI Ratio
vector<double> meanEIratioVector;
vector<double> meanInhEIratioVector;
vector<double> meanExcEIratioVector;
vector<double> standardDeviationEIratioVector;
vector<double> standardDeviationInhEIratioVector;
vector<double> standardDeviationExcEIratioVector;
vector<double> meanEIratioInfVector;
vector<double> meanInhEIratioInfVector;
vector<double> meanExcEIratioInfVector;
vector<double> standardDeviationEIratioInfVector;
vector<double> standardDeviationExcEIratioInfVector;
vector<double> standardDeviationInhEIratioInfVector;

//Threshold
vector<double> meanExcThresholdVector;
vector<double> meanInhThresholdVector;
vector<double> meanThresholdVector;
vector<double> sdExcThresholdVector;
vector<double> sdInhThresholdVector;
vector<double> sdThresholdVector;
vector<double> theta_exc_inf_vector;
vector<double> theta_inh_inf_vector;
vector<double> sdExcInfThresholdVector;
vector<double> sdInhInfThresholdVector;
vector<double> theta_inf_vector;
vector<double> sdInfThresholdVector;

//Mean Activities
vector<pair<double,double>> EM_dataVector_excitatory;
vector<pair<double,double>> EM_dataVector_inhibitory;
vector<double> gain_Exc_SD;
vector<double> gain_Inh_SD;
vector<double> m_exc_inf_vector;
vector<double> m_inh_inf_vector;
vector<double> m_exc_sd_inf_vector;
vector<double> m_inh_sd_inf_vector;

//total Inputs
vector<double> totalInputExcMeanVector;
vector<double> totalInputInhMeanVector;
vector<double> totalInputExcSDVector;
vector<double> totalInputInhSDVector;
vector<double> totalInputExcSDInfVector;
vector<double> totalInputInhSDInfVector;

//E Inputs
vector<double> EInputExcMeanVector;
vector<double> EInputInhMeanVector;
vector<double> EInputMeanVector;
vector<double> EInputExcMeanInfVector;
vector<double> EInputInhMeanInfVector;
vector<double> EInputMeanInfVector;
vector<double> EInputExcSDVector;
vector<double> EInputInhSDVector;
vector<double> EInputSDVector;
vector<double> EInputExcSDInfVector;
vector<double> EInputInhSDInfVector;
vector<double> EInputSDInfVector;

//I Inputs
vector<double> IInputExcMeanVector;
vector<double> IInputInhMeanVector;
vector<double> IInputMeanVector;
vector<double> IInputExcMeanInfVector;
vector<double> IInputInhMeanInfVector;
vector<double> IInputMeanInfVector;
vector<double> IInputExcSDVector;
vector<double> IInputInhSDVector;
vector<double> IInputSDVector;
vector<double> IInputExcSDInfVector;
vector<double> IInputInhSDInfVector;
vector<double> IInputSDInfVector;




//Simulation
Balanced_Network* network_scale;

for (int j=1;j<=Num_Scale1;j++){
for (int i=1;i<=Num_Scale2;i++){
  srand(6);
  network_scale = new
  Balanced_Network(Num_Excitatory_Neurons_Scale,Num_Inhibitory_Neurons_Scale
    , K_Scale,J_EE, J_EI, J_IE, J_II,m_0_Scale,
    externalRateFactor, j*0.1,i*0.02,update_times_Scale);
  for (int k=0;k<update_times_Scale;k++){
    Neuron* temp = network_scale->chooseRandomNeuron();
    network_scale->update2();
    //cout << 1;
  }

  //Phi & Lambda
  phiVector.push_back(network_scale->phi);
  lambdaVector.push_back(network_scale->lambda);

  //Mean Activity
  pair<double,double> dataPointExc = network_scale->getEM_data_exc2();
  pair<double,double> dataPointInh = network_scale->getEM_data_inh2();
  EM_dataVector_excitatory.push_back(dataPointExc);
  EM_dataVector_inhibitory.push_back(dataPointInh);
  gain_Exc_SD.push_back(network_scale->getEM_data_exc_sd());
  gain_Inh_SD.push_back(network_scale->getEM_data_inh_sd());
  m_exc_sd_inf_vector.push_back(network_scale->getM_exc_sd_inf());
  m_inh_sd_inf_vector.push_back(network_scale->getM_inh_sd_inf());
  m_exc_inf_vector.push_back(network_scale->getM_exc_inf());
  m_inh_inf_vector.push_back(network_scale->getM_inh_inf());

  //EI ratios
  network_scale->addEI_Ratios();
  meanEIratioVector.push_back(mean(network_scale->getEI_Ratios()));
  meanExcEIratioVector.push_back(mean(network_scale->getExcEI_Ratios()));
  meanInhEIratioVector.push_back(mean(network_scale->getInhEI_Ratios()));
  standardDeviationEIratioVector.push_back(
    standardDeviation(network_scale->getEI_Ratios()));
  standardDeviationExcEIratioVector.push_back(
    standardDeviation(network_scale->getExcEI_Ratios()));
  standardDeviationInhEIratioVector.push_back(standardDeviation(
    network_scale->getInhEI_Ratios()));
  meanEIratioInfVector.push_back(mean(network_scale->getEI_Ratios_inf()));
  meanExcEIratioInfVector.push_back(mean(network_scale->getExcEI_Ratios_inf()));
  meanInhEIratioInfVector.push_back(mean(network_scale->getInhEI_Ratios_inf()));
  standardDeviationEIratioInfVector.push_back(
      standardDeviation(network_scale->getEI_Ratios_inf()));
  standardDeviationExcEIratioInfVector.push_back(
      standardDeviation(network_scale->getExcEI_Ratios_inf()));
  standardDeviationInhEIratioInfVector.push_back(standardDeviation(
      network_scale->getInhEI_Ratios_inf()));

  //Threshold
  meanExcThresholdVector.push_back(network_scale->getMeanExcThreshold());
  meanInhThresholdVector.push_back(network_scale->getMeanInhThreshold());
  meanThresholdVector.push_back(network_scale->getMeanThreshold());
  sdExcThresholdVector.push_back(network_scale->getExcThresholdSD());
  sdInhThresholdVector.push_back(network_scale->getInhThresholdSD());
  sdThresholdVector.push_back(network_scale->getThresholdSD());
  theta_exc_inf_vector.push_back(network_scale->getTheta_exc_inf());
  theta_inh_inf_vector.push_back(network_scale->getTheta_inh_inf());
  theta_inf_vector.push_back(network_scale->getTheta_inf());
  sdExcInfThresholdVector.push_back(network_scale->getTheta_exc_inf_SD());
  sdInhInfThresholdVector.push_back(network_scale->getTheta_inh_inf_SD());
  sdInfThresholdVector.push_back(network_scale->getTheta_inf_SD());

  //Total Inputs
  totalInputExcMeanVector.push_back(network_scale->getTotalInputExcMean());
  totalInputInhMeanVector.push_back(network_scale->getTotalInputInhMean());
  totalInputExcSDVector.push_back(network_scale->getTotalInputExcSD());
  totalInputInhSDVector.push_back(network_scale->getTotalInputInhSD());
  totalInputExcSDInfVector.push_back(network_scale->getTotalInputExcSDInf());
  totalInputInhSDInfVector.push_back(network_scale->getTotalInputInhSDInf());

  //E Inputs
  EInputExcMeanVector.push_back(network_scale->getEInputExcMean());
  EInputInhMeanVector.push_back(network_scale->getEInputInhMean());
  EInputMeanVector.push_back(network_scale->getEInputMean());
  EInputExcMeanInfVector.push_back(network_scale->getEInputExcMeanInf());
  EInputInhMeanInfVector.push_back(network_scale->getEInputInhMeanInf());
  EInputMeanInfVector.push_back(network_scale->getEInputMeanInf());
  EInputExcSDVector.push_back(network_scale->getEInputExcSD());
  EInputInhSDVector.push_back(network_scale->getEInputInhSD());
  EInputSDVector.push_back(network_scale->getEInputSD());
  EInputExcSDInfVector.push_back(network_scale->getEInputExcSDInf());
  EInputInhSDInfVector.push_back(network_scale->getEInputInhSDInf());
  EInputSDInfVector.push_back(network_scale->getEInputSDInf());

  //I Inputs
  IInputExcMeanVector.push_back(network_scale->getIInputExcMean());
  IInputInhMeanVector.push_back(network_scale->getIInputInhMean());
  IInputMeanVector.push_back(network_scale->getIInputMean());
  IInputExcMeanInfVector.push_back(network_scale->getIInputExcMeanInf());
  IInputInhMeanInfVector.push_back(network_scale->getIInputInhMeanInf());
  IInputMeanInfVector.push_back(network_scale->getIInputMeanInf());
  IInputExcSDVector.push_back(network_scale->getIInputExcSD());
  IInputInhSDVector.push_back(network_scale->getIInputInhSD());
  IInputSDVector.push_back(network_scale->getIInputSD());
  IInputExcSDInfVector.push_back(network_scale->getIInputExcSDInf());
  IInputInhSDInfVector.push_back(network_scale->getIInputInhSDInf());
  IInputSDInfVector.push_back(network_scale->getIInputSDInf());

  //release memory
  delete network_scale;
}
}


//Outputting.

ofstream timeSeriesTxt1("Data/Time.txt");
ofstream timeSeriesTxt2("Data/TotalInput.txt");
ofstream timeSeriesTxt3("Data/ExcitatoryInput.txt");
ofstream timeSeriesTxt4("Data/InhibitoryInput.txt");

ofstream totalInputExcTxt("Data/TotalInputExcitatory.txt");
ofstream totalInputInhTxt("Data/TotalInputInhibitory.txt");
ofstream totalInputExcTimeTxt("Data/TotalInputExcitatoryTime.txt");
ofstream totalInputInhTimeTxt("Data/TotalInputInhibitoryTime.txt");

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
ofstream excEMSDTxt("Data/Exc_EM_SD.txt");
ofstream inhEMSDTxt("Data/Inh_EM_SD.txt");

ofstream thresholdTxt("Data/Threshold.txt");
ofstream spikeTimesTxt("Data/Spike_times.txt");

ofstream parametersTxt("Data/Parameters.txt");

ofstream excMeanThresholdTxt("Data/Exc_Mean_Threshold.txt");
ofstream excThresTimeTxt("Data/Exc_Thres_Time.txt");
ofstream inhMeanThresholdTxt("Data/Inh_Mean_Threshold.txt");
ofstream inhThresTimeTxt("Data/Inh_Thres_Time.txt");

ofstream adaptationGainTxt1("Data/Parameter_Scaling/Lambda.txt");
ofstream adaptationGainTxt2("Data/Parameter_Scaling/EI_Ratios_Mean.txt");
ofstream adaptationGainTxt3("Data/Parameter_Scaling/Phi.txt");
ofstream adaptationGainTxt4("Data/Parameter_Scaling/EI_Ratios_SD.txt");

ofstream plateauTimeTxt("Data/Plateau_Times.txt");
ofstream plateauNeuronsTxt("Data/Plateau_Neurons.txt");

ofstream adaptationGainTxt5("Data/Parameter_Scaling/Exc_MeanThreshold.txt");
ofstream adaptationGainTxt6("Data/Parameter_Scaling/Inh_MeanThreshold.txt");
ofstream adaptationGainTxt7("Data/Parameter_Scaling/MeanThreshold.txt");

ofstream adaptationGainTxt8("Data/Parameter_Scaling/Exc_ThresholdSD.txt");
ofstream adaptationGainTxt9("Data/Parameter_Scaling/Inh_ThresholdSD.txt");
ofstream adaptationGainTxt10("Data/Parameter_Scaling/ThresholdSD.txt");

ofstream adaptationGainTxt11("Data/Parameter_Scaling/Exc_EI_Ratio_Mean.txt");
ofstream adaptationGainTxt12("Data/Parameter_Scaling/Inh_EI_Ratio_Mean.txt");
ofstream adaptationGainTxt13("Data/Parameter_Scaling/Exc_EI_Ratio_SD.txt");
ofstream adaptationGainTxt14("Data/Parameter_Scaling/Inh_EI_Ratio_SD.txt");

ofstream adaptationGainTxt15("Data/Parameter_Scaling/Exc_M_inf.txt");
ofstream adaptationGainTxt16("Data/Parameter_Scaling/Inh_M_inf.txt");
ofstream adaptationGainTxt17("Data/Parameter_Scaling/Exc_Theta_inf.txt");
ofstream adaptationGainTxt18("Data/Parameter_Scaling/Inh_Theta_inf.txt");

ofstream adaptationGainTxt19("Data/Parameter_Scaling/TotalInput_Exc_Mean.txt");
ofstream adaptationGainTxt20("Data/Parameter_Scaling/TotalInput_Inh_Mean.txt");
ofstream adaptationGainTxt21("Data/Parameter_Scaling/TotalInput_Exc_SD.txt");
ofstream adaptationGainTxt22("Data/Parameter_Scaling/TotalInput_Inh_SD.txt");
ofstream adaptationGainTxt23("Data/Parameter_Scaling/TotalInput_Exc_SD_Inf.txt");
ofstream adaptationGainTxt24("Data/Parameter_Scaling/TotalInput_Inh_SD_Inf.txt");
ofstream adaptationGainTxt25("Data/Parameter_Scaling/Exc_ThresholdSD_inf.txt");
ofstream adaptationGainTxt26("Data/Parameter_Scaling/Inh_ThresholdSD_inf.txt");

ofstream adaptationGainTxt27("Data/Parameter_Scaling/Exc_EI_Ratio_Mean_inf.txt");
ofstream adaptationGainTxt28("Data/Parameter_Scaling/Inh_EI_Ratio_Mean_inf.txt");
ofstream adaptationGainTxt29("Data/Parameter_Scaling/EI_Ratio_Mean_inf.txt");
ofstream adaptationGainTxt30("Data/Parameter_Scaling/Exc_EI_Ratio_SD_inf.txt");
ofstream adaptationGainTxt31("Data/Parameter_Scaling/Inh_EI_Ratio_SD_inf.txt");
ofstream adaptationGainTxt32("Data/Parameter_Scaling/EI_Ratio_SD_inf.txt");

ofstream adaptationGainTxt33("Data/Parameter_Scaling/MeanThreshold_inf.txt");
ofstream adaptationGainTxt34("Data/Parameter_Scaling/ThresholdSD_inf.txt");

ofstream adaptationGainTxt35("Data/Parameter_Scaling/Exc_M_SD_inf.txt");
ofstream adaptationGainTxt36("Data/Parameter_Scaling/Inh_M_SD_inf.txt");

ofstream adaptationGainTxt37("Data/Parameter_Scaling/EInput_Exc_Mean.txt");
ofstream adaptationGainTxt38("Data/Parameter_Scaling/EInput_Inh_Mean.txt");
ofstream adaptationGainTxt39("Data/Parameter_Scaling/EInput_Mean.txt");
ofstream adaptationGainTxt40("Data/Parameter_Scaling/EInput_Exc_Mean_inf.txt");
ofstream adaptationGainTxt41("Data/Parameter_Scaling/EInput_Inh_Mean_inf.txt");
ofstream adaptationGainTxt42("Data/Parameter_Scaling/EInput_Mean_inf.txt");
ofstream adaptationGainTxt43("Data/Parameter_Scaling/EInput_Exc_SD.txt");
ofstream adaptationGainTxt44("Data/Parameter_Scaling/EInput_Inh_SD.txt");
ofstream adaptationGainTxt45("Data/Parameter_Scaling/EInput_SD.txt");
ofstream adaptationGainTxt46("Data/Parameter_Scaling/EInput_Exc_SD_inf.txt");
ofstream adaptationGainTxt47("Data/Parameter_Scaling/EInput_Inh_SD_Inf.txt");
ofstream adaptationGainTxt48("Data/Parameter_Scaling/EInput_SD_inf.txt");

ofstream adaptationGainTxt49("Data/Parameter_Scaling/IInput_Exc_Mean.txt");
ofstream adaptationGainTxt50("Data/Parameter_Scaling/IInput_Inh_Mean.txt");
ofstream adaptationGainTxt51("Data/Parameter_Scaling/IInput_Mean.txt");
ofstream adaptationGainTxt52("Data/Parameter_Scaling/IInput_Exc_Mean_inf.txt");
ofstream adaptationGainTxt53("Data/Parameter_Scaling/IInput_Inh_Mean_inf.txt");
ofstream adaptationGainTxt54("Data/Parameter_Scaling/IInput_Mean_inf.txt");
ofstream adaptationGainTxt55("Data/Parameter_Scaling/IInput_Exc_SD.txt");
ofstream adaptationGainTxt56("Data/Parameter_Scaling/IInput_Inh_SD.txt");
ofstream adaptationGainTxt57("Data/Parameter_Scaling/IInput_SD.txt");
ofstream adaptationGainTxt58("Data/Parameter_Scaling/IInput_Exc_SD_inf.txt");
ofstream adaptationGainTxt59("Data/Parameter_Scaling/IInput_Inh_SD_Inf.txt");
ofstream adaptationGainTxt60("Data/Parameter_Scaling/IInput_SD_inf.txt");


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

//Total Inputs network
for (int i=0;i<totalInputExc.size();i++){
  totalInputExcTimeTxt << totalInputExc[i].first << endl;
  totalInputExcTxt << totalInputExc[i].second << endl;
}

for (int i=0;i<totalInputInh.size();i++){
  totalInputInhTimeTxt << totalInputInh[i].first <<endl;
  totalInputInhTxt << totalInputInh[i].second << endl;
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
for (int i=0;i<exc_mean.size();i++){
  excTimeTxt << exc_mean[i].first << endl;
  excMean_ActivityTxt << (exc_mean[i].second/Num_Excitatory_Neurons) << endl;
}

for (int i=0;i<inh_mean.size();i++){
  inhTimeTxt << inh_mean[i].first << endl;
  inhMean_ActivityTxt << (inh_mean[i].second/Num_Inhibitory_Neurons) << endl;
}


//External Input vs. Mean activity
for (int i=0;i<EM_dataVector_excitatory.size();i++){
  excEMRateTxt << EM_dataVector_excitatory[i].first << endl;
  excEMAtvTxt << EM_dataVector_excitatory[i].second << endl;
  inhEMRateTxt << EM_dataVector_inhibitory[i].first << endl;
  inhEMAtvtxt << EM_dataVector_inhibitory[i].second << endl;
  excEMSDTxt << gain_Exc_SD[i] << endl;
  inhEMSDTxt << gain_Inh_SD[i]<< endl;
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


//Parameter Scaling

for (int i=0;i<phiVector.size();i++){
  adaptationGainTxt1 << lambdaVector[i] << endl;
  adaptationGainTxt2 << meanEIratioVector[i] << endl;
  adaptationGainTxt3 << phiVector[i]<< endl;
  adaptationGainTxt4 << standardDeviationEIratioVector[i] << endl;
  adaptationGainTxt5 << meanExcThresholdVector[i] << endl;
  adaptationGainTxt6 << meanInhThresholdVector[i] << endl;
  adaptationGainTxt7 << meanThresholdVector[i]<< endl;
  adaptationGainTxt8 << sdExcThresholdVector[i] << endl;
  adaptationGainTxt9 << sdInhThresholdVector[i] << endl;
  adaptationGainTxt10 << sdThresholdVector[i] << endl;
  adaptationGainTxt11 << meanExcEIratioVector[i]<<endl;
  adaptationGainTxt12 << meanInhEIratioVector[i] << endl;
  adaptationGainTxt13 << standardDeviationExcEIratioVector[i] << endl;
  adaptationGainTxt14 << standardDeviationInhEIratioVector[i] << endl;
  adaptationGainTxt15 << m_exc_inf_vector[i]<< endl;
  adaptationGainTxt16 << m_inh_inf_vector[i]<< endl;
  adaptationGainTxt17 << theta_exc_inf_vector[i] << endl;
  adaptationGainTxt18 << theta_inh_inf_vector[i] << endl;
  adaptationGainTxt19 << totalInputExcMeanVector[i] << endl;
  adaptationGainTxt20 << totalInputInhMeanVector[i] << endl;
  adaptationGainTxt21 << totalInputExcSDVector[i] << endl;
  adaptationGainTxt22 << totalInputInhSDVector[i] << endl;
  adaptationGainTxt23 << totalInputExcSDInfVector[i] << endl;
  adaptationGainTxt24 << totalInputInhSDInfVector[i] << endl;
  adaptationGainTxt25 << sdExcInfThresholdVector[i] << endl;
  adaptationGainTxt26 << sdInhInfThresholdVector[i] << endl;

  adaptationGainTxt27 << meanEIratioInfVector[i] << endl;
  adaptationGainTxt28 << meanExcEIratioInfVector[i] << endl;
  adaptationGainTxt29 << meanInhEIratioInfVector[i] << endl;
  adaptationGainTxt30 << standardDeviationEIratioInfVector[i]<<endl;
  adaptationGainTxt31 << standardDeviationExcEIratioInfVector[i] <<endl;
  adaptationGainTxt32 << standardDeviationInhEIratioInfVector[i] << endl;

  adaptationGainTxt33 << theta_inf_vector[i] <<endl;
  adaptationGainTxt34 << sdInfThresholdVector[i] << endl;

  adaptationGainTxt35 << m_exc_sd_inf_vector[i] << endl;
  adaptationGainTxt36 << m_inh_sd_inf_vector[i] <<endl;

  adaptationGainTxt37 << EInputExcMeanVector[i] <<endl;
  adaptationGainTxt38 << EInputInhMeanVector[i] <<endl;
  adaptationGainTxt39 << EInputMeanVector[i] <<endl;
  adaptationGainTxt40 << EInputExcMeanInfVector[i] <<endl;
  adaptationGainTxt41 << EInputInhMeanInfVector[i] <<endl;
  adaptationGainTxt42 << EInputMeanInfVector[i] <<endl;
  adaptationGainTxt43 << EInputExcSDVector[i] <<endl;
  adaptationGainTxt44 << EInputInhSDVector[i] <<endl;
  adaptationGainTxt45 << EInputSDVector[i] <<endl;
  adaptationGainTxt46 << EInputExcSDInfVector[i] <<endl;
  adaptationGainTxt47 << EInputInhSDInfVector[i] <<endl;
  adaptationGainTxt48 << EInputSDInfVector[i] <<endl;

  adaptationGainTxt49 << IInputExcMeanVector[i] <<endl;
  adaptationGainTxt50 << IInputInhMeanVector[i] <<endl;
  adaptationGainTxt51 << IInputMeanVector[i] <<endl;
  adaptationGainTxt52 << IInputExcMeanInfVector[i] <<endl;
  adaptationGainTxt53 << IInputInhMeanInfVector[i] <<endl;
  adaptationGainTxt54 << IInputMeanInfVector[i] <<endl;
  adaptationGainTxt55 << IInputExcSDVector[i] <<endl;
  adaptationGainTxt56 << IInputInhSDVector[i] <<endl;
  adaptationGainTxt57 << IInputSDVector[i] <<endl;
  adaptationGainTxt58 << IInputExcSDInfVector[i] <<endl;
  adaptationGainTxt59 << IInputInhSDInfVector[i] <<endl;
  adaptationGainTxt60 << IInputSDInfVector[i] <<endl;
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
parametersTxt << "External Input: " + to_string(nVector[0]->externalInput) << endl;
parametersTxt << "Number of updates: " + to_string(update_times) <<endl;
parametersTxt << "Adaptation Jump: " + to_string(nVector[0]->adaptation_jump) << endl;
parametersTxt << "lambda: " + to_string(nVector[0]->decay_constant) << endl;
parametersTxt << " " << endl;
parametersTxt << "Scaling part: " << endl;
parametersTxt << "Number of networks: " + to_string(Num_Scale2*Num_Scale1) << endl;
parametersTxt << "Number of neurons: " + to_string(Num_Excitatory_Neurons_Scale+ Num_Inhibitory_Neurons_Scale) << endl;
parametersTxt << "K: " + to_string(K_Scale) << endl;
parametersTxt << "Number of updates: " + to_string(update_times_Scale) << endl;
parametersTxt << "m_0_Scale: " + to_string(m_0_Scale) << endl;

delete neural_network;

  return 0;
}
