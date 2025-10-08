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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "RunAction.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "SteppingActionMessenger.hh"
#include "G4AutoLock.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include <chrono>

SteppingAction::SteppingAction(EventAction* eventAction, RunAction* RuAct)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fRunAction(RuAct),
  fBackscatterFilename(),
  fSteppingMessenger(),
  fCollectionAltitude(450.0)
{
  fSteppingMessenger = new SteppingActionMessenger(this);
}

SteppingAction::~SteppingAction(){delete fSteppingMessenger;}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Dividing by a unit outputs data in that unit, so divisions by keV result in outputs in keV
  // https://geant4-internal.web.cern.ch/sites/default/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
  G4Track* track = step->GetTrack();
  const G4int trackID = track->GetTrackID();
  const G4ThreeVector position = track->GetPosition();
  const G4ThreeVector momentumDirection = track->GetMomentumDirection();
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  const G4double preStepKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  const G4double postStepKineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  const G4double trackWeight = track->GetWeight();
  
  // ===========================
  // Guard Block
  // ===========================

  // Check for NaN energy
  if(std::isnan(postStepKineticEnergy))
  {  
    G4cout << "WARNING: Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: NaN energy. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }
  // Check for exceeding 1 second of simulation time
  if(track->GetProperTime()/second > 1)
  {
    G4cout << "WARNING: Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: Exceeded 1s simulation time. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }
  // Check for stuck photons. Occassionally they seem to get 'wedged' between atmospheric layers and stop propagating without being automatically killed, hanging the program forever
  if((step->GetStepLength()/m < 1e-12) && particleName == "gamma"){
    G4cout << "WARNING: Killed " << particleName << " at " << postStepKineticEnergy/keV << " keV. Reason: Stuck gamma. Current step length: " << step->GetStepLength()/m << " m" << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // ===========================
  // Begin Data Logging
  // ===========================
  // Energy Spectrum Tracking






  if(particleName != "gamma"){track->SetTrackStatus(fStopAndKill);}
  







  // Get altitude information
  G4double preStepAlt_km  = (step->GetPreStepPoint()->GetPosition().z()/km) + 500.0;
  G4double postStepAlt_km = (step->GetPostStepPoint()->GetPosition().z()/km) + 500.0;
  
  // Get altitude indices (float) of start and stop point
  G4double preStepAltitudeIndex = (preStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;
  G4double postStepAltitudeIndex = (postStepAlt_km -  fRunAction->fMinSampleAltitude_km) / fRunAction->altitudeSpacing_km;

  // Kick out if step is entirely outside the altitudes we care about
  bool preStepInRange = (0 <= preStepAltitudeIndex) && (preStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  bool postStepInRange = (0 <= postStepAltitudeIndex) && (postStepAltitudeIndex < (fRunAction->fNumberOfSamplePlanes-1));
  if( (preStepInRange == false) && (postStepInRange == false) ){ return; }
  
  // Get bounding indices of planes that have been crossed
  int startIdx = std::ceil(std::min(preStepAltitudeIndex, postStepAltitudeIndex));
  int stopIdx = std::floor(std::max(preStepAltitudeIndex, postStepAltitudeIndex));










  G4cout.precision(15);
  G4double dz = std::abs(postStepAlt_km - preStepAlt_km);
  if(dz < 0.5){
    G4cout << std::fixed << "(" << trackID << ") " << preStepAlt_km << " -> " << postStepAlt_km << ", dz = " << dz << ", idxs = " << startIdx << " -> " << stopIdx << G4endl;
  }


  

  





  // Kick out particles that didn't cross any planes
  if(startIdx > stopIdx){return;}
  
  // Loop over crossed planes and add to energy spectra
  for(int altitudeIndex = startIdx; altitudeIndex <= stopIdx; altitudeIndex++){
    // Kick out invalid indices
    if( (altitudeIndex < 0) || (altitudeIndex > (fRunAction->fNumberOfSamplePlanes-1)) ){continue;}
    G4cout << "logging " << altitudeIndex << G4endl;



    // Do an interpolation to get approximate energy at the plane crossing
    G4double t = (fRunAction->sampleAltitudes_km[altitudeIndex] - preStepAlt_km) / (postStepAlt_km - preStepAlt_km);
    G4double crossingEnergy = preStepKineticEnergy + (t * (postStepKineticEnergy - preStepKineticEnergy)); // Linear interpolation

    // Find the energy bin this particle resides in utilizing regular spacing to directly calculate the index.
    int energyIndex = std::floor(logbase(fRunAction->histogramFactor, (crossingEnergy/keV)/(fRunAction->fEnergyMinkeV)));
    if( (energyIndex < 0) || (energyIndex > (fRunAction->fNumberOfEnergyBins-1)) ){
      continue; // Don't record energy if particle energy is out of range of histogram
    }
    
    // Add count to histogram based on particle type
    if(particleName == "e-")    {fRunAction->electronCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "proton"){fRunAction->protonCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "gamma") {fRunAction->gammaCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "alpha") {fRunAction->alphaCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
  }
}

G4double SteppingAction::logbase(G4double base, G4double x){
  // Log base of x, log_base(x)
  return log2(x)/log2(base);
}