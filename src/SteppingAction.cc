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





#include <iostream>
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
  
  // Get altitude information
  G4double preStepAlt_km  = (step->GetPreStepPoint()->GetPosition().z()/km) + 500.0;
  G4double postStepAlt_km = (step->GetPostStepPoint()->GetPosition().z()/km) + 500.0;
  
  // Kick out particles outside the range of altitudes we care about
  bool preStepInRange = (preStepAlt_km <= fRunAction->fMaxSampleAltitude_km) && (preStepAlt_km >= fRunAction->fMinSampleAltitude_km);
  bool postStepInRange = (postStepAlt_km <= fRunAction->fMaxSampleAltitude_km) && (postStepAlt_km >= fRunAction->fMinSampleAltitude_km);
  if( (preStepInRange == false) && (postStepInRange == false) ){ return; }




  //int n_planes_crossed = 0;




  // Loop over every sample altitude and see if this particle has crossed that plane
  for(int altitudeIndex = 0; altitudeIndex <  fRunAction->fNumberOfSamplePlanes; altitudeIndex++){
    bool crossedPlane = (preStepAlt_km > fRunAction->sampleAltitudes_km[altitudeIndex]) != (postStepAlt_km > fRunAction->sampleAltitudes_km[altitudeIndex]);
    if( crossedPlane == false ){continue;} // Don't proceed if we haven't crossed the plane
  
    /*
    n_planes_crossed++;
    if(n_planes_crossed > 1){
      G4cout << "Crossed " << n_planes_crossed << " planes." << G4endl;
    }
    */

    // Get energy at plane crossing point
    G4double crossingEnergy = postStepKineticEnergy;

    // If it's not safe to use pre- and post- kinetic energy interchangably, do an interpolation to get approximate energy at the plane crossing.
    if( std::abs((postStepKineticEnergy - preStepKineticEnergy)/preStepKineticEnergy) > 1e-3 ){
      G4double t = (fRunAction->sampleAltitudes_km[altitudeIndex] - preStepAlt_km) / (postStepAlt_km - preStepAlt_km);
      crossingEnergy = preStepKineticEnergy + (t * (postStepKineticEnergy - preStepKineticEnergy)); // Linear interpolation
    }

    // Find the energy bin this particle resides in. Do NOT do a comparison to every bin - utilize the fact that bins are logspaced
    // to directly calculate the index.
    int energyIndex = std::floor(logbase(fRunAction->histogramFactor, (crossingEnergy/keV)/(fRunAction->fEnergyMinkeV)));
    if( (energyIndex < 0) || (energyIndex > (fRunAction->fNumberOfEnergyBins-1)) ){
      continue; // Don't record energy if particle energy is out of range of histogram
    }
    
    // Add count to histogram based on particle type
    if(particleName == "e-")    {fRunAction->electronCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "proton"){fRunAction->protonCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "gamma") {fRunAction->gammaCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    if(particleName == "alpha") {fRunAction->alphaCounts[altitudeIndex][energyIndex] += 1 * trackWeight;}
    
    // If the particle's altitude change is less than half the space between altitude planes, it is not possible to have crossed
    // multiple planes in this step and we can thus break the loop for speed.
    if( std::abs(postStepAlt_km - postStepAlt_km) < fRunAction->altitudeSpacing_km/2){
      break;
    }
  }
}

G4double SteppingAction::logbase(G4double base, G4double x){
  // Log base of x, log_base(x)
  return log2(x)/log2(base);
}