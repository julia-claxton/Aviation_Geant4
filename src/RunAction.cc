//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4Transportation.hh"
#include "G4CoupledTransportation.hh"
#include "G4Electron.hh"

#include "myHistogram.hh"
#include "RunActionMessenger.hh"

#include <fstream>
#include <regex>
#include <filesystem>
#include "../include/csv.h" // For quickly parsing .csv files. Author credit in file.

RunAction::RunAction():
  G4UserRunAction(),
  fRunActionMessenger(),
  fEnergySpectraFileName()
  
{
  // Set killing energies
  fWarningEnergy = 0.01 * keV; // Particles below this energy are killed after 1 step. Value arbitrary 
  fImportantEnergy = 0.1 * keV; // Particles above this energy are killed after fNumberOfTrials if they are looping. Value arbitrary 
  fNumberOfTrials = 1000; // Number of trials before a looping 'important' particle is killed. Value arbitrary

  fRunActionMessenger = new RunActionMessenger(this); 

  // Create sample planes
  sampleAltitudes_km = linspace(fMinSampleAltitude_km, fMaxSampleAltitude_km, fNumberOfSamplePlanes);

  // Create energy bins
  energyBinEdges_keV = linspace(std::log10(fEnergyMinkeV), std::log10(fEnergyMaxkeV), fNumberOfEnergyBins + 1); 
  for(int i = 0; i < fNumberOfEnergyBins+1; i++){
    energyBinEdges_keV.at(i) = pow(10, energyBinEdges_keV.at(i));
  }
  
  // Create histograms initialized to zero
  protonCounts.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));
  electronCounts.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));
  gammaCounts.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));
  alphaCounts.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));
  
  // Precalculate things
  histogramFactor = pow(10, (std::log10(fEnergyMaxkeV) - std::log10(fEnergyMinkeV)) / fNumberOfEnergyBins);
  altitudeSpacing_km = std::abs(sampleAltitudes_km[1] - sampleAltitudes_km[0]);
}

RunAction::~RunAction()
{
  delete fRunActionMessenger;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // If we are the main thread, create backscatter file and write header
  int threadID = G4Threading::G4GetThreadId();
  if(threadID == -1)
  {
    // First, make sure that the build directory set by the user is correct
    std::filesystem::path resultsPath = fEnergySpectraFileName.c_str();
    std::filesystem::path buildDirectory = resultsPath.parent_path().parent_path();
    if(std::filesystem::is_directory(buildDirectory) == false)
    {
      G4cerr << G4endl << "*** ERROR: User-specified build directory " << buildDirectory << " does not exist. This path is user-specified in set_simulation_parameters.mac. Check that G4EPP_BUILD_DIR in set_simulation_parameters.mac matches your build directory and does not have a slash at the end." << G4endl << G4endl;
      throw;
    }
  }
  // Otherwise, print startup message
  else
  {
    // Pad with spaces to have consistent print location
    int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
    int paddingLength = std::to_string(nThreads).length() - std::to_string(threadID).length();

    // Get current time to print
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);

    // Print startup message
    G4cout << std::string(paddingLength, ' ') << "(" << std::put_time(&tm, "%F %T") <<") STARTING: Thread " << threadID << G4endl;
  }

  // Change parameters for looping particles
  ChangeLooperParameters( G4Electron::Definition() );
}

void RunAction::ChangeLooperParameters(const G4ParticleDefinition* particleDef)
{
  if(particleDef == nullptr)
    particleDef = G4Electron::Definition();
  auto transportPair= findTransportation(particleDef);
  auto transport = transportPair.first;
  auto coupledTransport = transportPair.second;

  if(transport != nullptr)
  { 
    // Change the values of the looping particle parameters of Transportation 
    if(fWarningEnergy >= 0.0)
      transport->SetThresholdWarningEnergy(  fWarningEnergy ); 
    if(fImportantEnergy >= 0.0)
      transport->SetThresholdImportantEnergy(  fImportantEnergy ); 
    if(fNumberOfTrials > 0)
      transport->SetThresholdTrials( fNumberOfTrials );
  }
  else if(coupledTransport != nullptr)
  { 
    // Change the values for Coupled Transport
    if(fWarningEnergy >= 0.0)
      coupledTransport->SetThresholdWarningEnergy(fWarningEnergy); 
    if(fImportantEnergy >= 0.0)
      coupledTransport->SetThresholdImportantEnergy(fImportantEnergy); 
    if(fNumberOfTrials > 0)
      coupledTransport->SetThresholdTrials(fNumberOfTrials);
  }
}

void RunAction::EndOfRunAction(const G4Run*)
{
  // Get thread ID to see if we are main thread or not
  int threadID = G4Threading::G4GetThreadId();

  // If we are not the main thread, write energy deposition and backscatter to file and exit
  if(threadID != -1)
  {
    std::string particlesToWrite[] = {"electron", "proton", "gamma", "alpha"};
    std::vector<std::vector<std::vector<G4double>>> dataToWrite = {electronCounts, protonCounts, gammaCounts, alphaCounts};

    for(int particleIndex = 0; particleIndex < 4; particleIndex++){
      // Write energy deposition to file
      std::string filepath = 
        fEnergySpectraFileName.substr(0, fEnergySpectraFileName.length()-4) 
        + "_" + particlesToWrite[particleIndex]
        + "_thread" + std::to_string(threadID)
      + ".csv";

      writeThreadHistogramToFile(filepath, dataToWrite[particleIndex]);
    }

    // Done writing data, now print status message
    // Pad with spaces to have consistent print location
    int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
    int paddingLength = std::to_string(nThreads).length() - std::to_string(threadID).length();
    
    // Get current time to print
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);

    // Print message
    G4cout << std::string(paddingLength, ' ') << "(" << std::put_time(&tm, "%F %T") <<") \033[0;32mFINISHED: Thread " << threadID << "\033[0m" << G4endl;
    return;
  }

  // If we are the main thread, merge energy deposition datafiles from each thread. Main thread ends after workers are done, so this is the end of the simulation
  G4cout << "Merging thread-specific data... ";

  std::string particlesToWrite[] = {"electron", "proton", "gamma", "alpha"};
  





  // TODO loop this over thread and particle type
  for(int particleIndex = 0; particleIndex < 4; particleIndex++){
    mainFilename = fEnergySpectraFileName.substr(0, fEnergySpectraFileName.length()-4) + particlesToWrite[particleIndex] + ".csv";

    int thread = 0;

    std::string threadFilepath = 
      fEnergySpectraFileName.substr(0, fEnergySpectraFileName.length()-4) 
      + "_" + particlesToWrite[particleIndex]
      + "_thread" + std::to_string(thread)
    + ".csv";

    // csv.h seems to be hard to use when you have hundreds of columns so I have to resort to writing my own reader. Lays down facedown on the floor.
    std::vector<std::vector<G4double>> threadData = readThreadFile(threadFilepath);
    
    // TODO delete threadfile
    
    // end loop
    G4cout << "TODO" << G4endl; throw;
    
    writeMainHistogramToFile(fEnergySpectraFileName, threadData); // TEMP for testing
  }

  G4cout << "Done" << G4endl;
}

std::vector<G4double> RunAction::linspace(G4double start, G4double stop, int n){
  std::vector<G4double> result;
  G4double stepSize = (stop - start) / (n-1);

  for(int i = 0; i < n; i++){
    result.push_back(start + (i * stepSize));
  }
  return result;
}

void RunAction::writeMainHistogramToFile(std::string filename, std::vector<std::vector<G4double>> histogram)
{
  // Open file
  std::ofstream dataFile;
  dataFile.open(filename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Write header
  dataFile << "altitude_km,";
  for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
    dataFile << energyBinEdges_keV[energyIndex] << "keV-" << energyBinEdges_keV[energyIndex+1] << "keV";
    energyIndex == fNumberOfEnergyBins-1 ? dataFile << "\n" : dataFile << ",";
  }

  // Write rows
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
    // Write altitude label
    dataFile << sampleAltitudes_km[altitudeIndex] << ",";

    // Write energy spectrum
    for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
      dataFile << histogram[altitudeIndex][energyIndex];
      energyIndex == fNumberOfEnergyBins-1 ? dataFile << "\n" : dataFile << ",";
    }
  }

  // Close file
  dataFile.close();
}

void RunAction::writeThreadHistogramToFile(std::string filename, std::vector<std::vector<G4double>> histogram)
{
  // Open file
  std::ofstream dataFile;
  dataFile.open(filename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // No header, go right to the rows
  for(int altitudeIndex = 0; altitudeIndex < fNumberOfSamplePlanes; altitudeIndex++){
    for(int energyIndex = 0; energyIndex < fNumberOfEnergyBins; energyIndex++){
      // Write energy spectrum
      dataFile << histogram[altitudeIndex][energyIndex];
      energyIndex == fNumberOfEnergyBins-1 ? dataFile << "\n" : dataFile << ",";
    }
  }

  // Close file
  dataFile.close();
}

std::vector<std::vector<G4double>> RunAction::readThreadFile(std::string filename){
  // Allocate result
  std::vector<std::vector<G4double>> result;
  result.resize(fNumberOfSamplePlanes, std::vector<G4double>(fNumberOfEnergyBins, 0));

  // Open file
  std::ifstream file;
  file.open(filename, std::ios_base::out); // Open file in write mode to overwrite any previous results

  // Parse lines
  int dim1Index = 0;
  int dim2Index = 0;

  std::string line;
  std::string token;
  while( std::getline(file, line) ){
    std::istringstream word(line);
    while ( std::getline(word, token, ',') ){
      result[dim1Index][dim2Index] = std::stod(token);
      dim2Index++;
    }
    dim1Index++;
    dim2Index = 0;
  }
  file.close();
  return result;
}

std::pair<G4Transportation*, G4CoupledTransportation*> RunAction::findTransportation(const G4ParticleDefinition* particleDef, bool reportError)
{
  const auto *partPM = particleDef->GetProcessManager();
    
  G4VProcess* partTransport = partPM->GetProcess("Transportation");
  auto transport= dynamic_cast<G4Transportation*>(partTransport);

  partTransport = partPM->GetProcess("CoupledTransportation");
  auto coupledTransport = dynamic_cast<G4CoupledTransportation*>(partTransport);

  if(reportError && !transport && !coupledTransport)
  {
    G4cerr << "Unable to find Transportation process for particle type "
           << particleDef->GetParticleName()
           << "  ( PDG code = " << particleDef->GetPDGEncoding() << " ) "
    << G4endl;
  }
  
  return std::make_pair( transport, coupledTransport );
}
