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
  fEnergyDepositionFileName()
  
{
  fWarningEnergy = 0.01 * keV; // Particles below this energy are killed after 1 step. Value arbitrary 
  fImportantEnergy = 0.1 * keV; // Particles above this energy are killed after fNumberOfTrials if they are looping. Value arbitrary 
  fNumberOfTrials = 1000; // Number of trials before a looping 'important' particle is killed. Value arbitrary

  fRunActionMessenger = new RunActionMessenger(this); 

  // Set up energy spectrum histograms
  G4cout << fNumberEnergyBins << " bins in spectrum" << G4endl;

  G4double energyMinMeV = 1e-3; // TEMP
  G4double energyMaxMeV = 1e3; // TEMP

  // TODO logspace
  std::vector<G4double> energyBinEdges = linspace(energyMinMeV, energyMaxMeV, fNumberEnergyBins + 1); 

  // todo move fNumberEnergyBins inside here, or move min/max out to class def in .hh
  // todo construct the full 2d histogram (x3 for each species w/ stopping power data)
  // todo get those full histograms into global runaction scope so steppingaction can add to it

  throw;
  
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
    std::filesystem::path resultsPath = fEnergyDepositionFileName.c_str();
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
  /*
  // Get thread ID to see if we are main thread or not
  int threadID = G4Threading::G4GetThreadId();

  // If we are not the main thread, write energy deposition and backscatter to file and exit
  if(threadID != -1)
  {
    // Write energy deposition to file
    std::string energyDepositionThreadFilename = fEnergyDepositionFileName.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadID) + ".csv"; // Thread-specific filename
    fEnergyDepositionHistogram->WriteHistogramToFile(energyDepositionThreadFilename);

    // Write backscatter to file
    std::string backscatterThreadFilename = std::regex_replace(fEnergyDepositionFileName, std::regex("energy_deposition"), "backscatter");
    backscatterThreadFilename = backscatterThreadFilename.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadID) + ".csv"; // Thread-specific filename
    writeBackscatterToFile(backscatterThreadFilename);

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

  // Create main energy deposition data structure
  myHistogram* mainEnergyDepositionHistogram = new myHistogram(); // 1000 km in 1 km bins

  // Create main backscatter file
  std::ofstream backscatterFile;
  std::string backscatterFilename = std::regex_replace(fEnergyDepositionFileName, std::regex("energy_deposition"), "backscatter");
  backscatterFile.open(backscatterFilename, std::ios_base::out); // Open file in write mode to overwrite any previous results
  backscatterFile << "particle_name,particle_weight,particle_energy_keV,particle_pitch_angle_deg,momentum_direction_x,momentum_direction_y,momentum_direction_z,x_meters,y_meters,z_meters\n";

  int nThreads = G4Threading::GetNumberOfRunningWorkerThreads();
  for(int threadFileToMerge = 0; threadFileToMerge < nThreads; threadFileToMerge++)
  {
    // ENERGY DEPOSITION:
    std::string energyDepositionThreadFilename = fEnergyDepositionFileName.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadFileToMerge) + ".csv"; // Thread-specific filename

    // Read in energy deposition from this thread via csv
    io::CSVReader<2> in(energyDepositionThreadFilename);
    in.read_header(io::ignore_extra_column, "altitude_km", "energy_deposition_kev");
    int altitudeAddress; double energy_deposition;

    // For each row in the file, add energy deposition to main histogram
    while(in.read_row(altitudeAddress, energy_deposition)){ 
      mainEnergyDepositionHistogram->AddCountToBin(altitudeAddress, energy_deposition);
    }
    // Delete this thread-specific file
    std::remove(energyDepositionThreadFilename.c_str());


    // BACKSCATTER:
    std::string backscatterThreadFilename = std::regex_replace(fEnergyDepositionFileName, std::regex("energy_deposition"), "backscatter");
    backscatterThreadFilename = backscatterThreadFilename.substr(0, fEnergyDepositionFileName.length()-4) + "_thread" + std::to_string(threadFileToMerge) + ".csv"; // Thread-specific filename

    // Move on if file doesn't exist (i.e. no backscatter from this thread)
    if(std::filesystem::exists(backscatterThreadFilename) == false){continue;}

    // Read in backscatter from this thread via csv
    io::CSVReader<10> inBackscatter(backscatterThreadFilename);
    inBackscatter.read_header(io::ignore_extra_column, "particle_name", "particle_weight", "particle_energy_keV", "particle_pitch_angle_deg", "momentum_direction_x", "momentum_direction_y", "momentum_direction_z", "x_meters", "y_meters", "z_meters");
    std::string particle_name;
    double particle_weight;
    double particle_energy_keV;
    double particle_pitch_angle_deg;
    double momentum_direction_x;
    double momentum_direction_y;
    double momentum_direction_z;
    double x_meters;
    double y_meters;
    double z_meters;

    // For each row in the file, add energy deposition to main histogram
    while(inBackscatter.read_row(particle_name, particle_weight, particle_energy_keV, particle_pitch_angle_deg, momentum_direction_x, momentum_direction_y, momentum_direction_z, x_meters, y_meters, z_meters)){ 
      backscatterFile << 
        particle_name            << "," <<
        particle_weight          << "," <<
        particle_energy_keV      << "," <<
        particle_pitch_angle_deg << "," <<
        momentum_direction_x     << "," <<
        momentum_direction_y     << "," <<
        momentum_direction_z     << "," <<
        x_meters                 << "," <<
        y_meters                 << "," <<
        z_meters                 << "\n"
      ;
    }

    // Delete the thread-specific file
    std::remove(backscatterThreadFilename.c_str());
  }
  mainEnergyDepositionHistogram->WriteHistogramToFile(fEnergyDepositionFileName);
  backscatterFile.close();
  */

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
