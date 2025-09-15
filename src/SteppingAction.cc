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

// Initialize autolock for multiple threads writing into a single file
namespace{ G4Mutex aMutex=G4MUTEX_INITIALIZER; } 

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
  // ===========================
  // Guard Block
  // ===========================
  G4Track* track = step->GetTrack();
  const G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4double preStepKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  G4double postStepKineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  G4double trackWeight = track->GetWeight();

  // Check for NaN energy
  if(std::isnan(postStepKineticEnergy))
  {  
    G4cout << "Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: NaN energy. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // Check for exceeding 1 second of simulation time
  if(track->GetProperTime()/second > 1)
  {
    G4cout << "Killed " << particleName << "at " << postStepKineticEnergy/keV << " keV. Reason: Exceeded 1s simulation time. Process: " << step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // Check for stuck photons. Occassionally they seem to get 'wedged' between atmospheric layers and stop propagating without being automatically killed, hanging the program forever
  if((step->GetStepLength()/m < 1e-12) && particleName == "gamma"){
    G4cout << "Killed " << particleName << " at " << postStepKineticEnergy/keV << " keV. Reason: Stuck gamma. Current step length: " << step->GetStepLength()/m << " m" << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  // ===========================
  // Begin Data Logging
  // ===========================
  // Dividing by a unit outputs data in that unit, so divisions by keV result in outputs in keV
  // https://geant4-internal.web.cern.ch/sites/default/files/geant4/collaboration/working_groups/electromagnetic/gallery/units/SystemOfUnits.html
  const G4ThreeVector position = track->GetPosition();
  const G4ThreeVector momentumDirection = track->GetMomentumDirection();

  // ===========================
  // Energy Deposition Tracking
  // ===========================
  // Add energy deposition to histogram owned by RunAction, which is written to a results file at the end of each thread's simulation run  
  G4double zPos = position.z(); // Particle altitude in world coordinates
  
  // Kill particles that go below ground level
  if((500.0 + (zPos/km)) < 0){ track->SetTrackStatus(fStopAndKill); }

  // Log energy deposition
  G4int altitudeAddress = std::floor(500.0 + zPos/km); // Index to write to. Equal to altitude above sea level in km, to lowest whole km
  if(altitudeAddress > 0 && altitudeAddress < 1000) 
  {
    G4double weightedEnergyDeposition = step->GetTotalEnergyDeposit() * trackWeight;
    LogEnergy(altitudeAddress, weightedEnergyDeposition/keV); // Threadlocking occurs inside LogEnergy
  }

  // ===========================
  // Backscatter Tracking
  // ===========================
  // Write backscatter to memory if detected
  // Backscatter is defined as a particle above the collection altitude and moving upwards in world coordinates
  // We subtract 500 from the collection altitude because +500.0 km above sea level ==> z = 0.0 in world coordinates
  if( (position.z()/km > fCollectionAltitude-500.0) && (momentumDirection.z() > 0) ) 
  {
    // Calculate particle pitch angle
    G4double spacetimePoint[4] = {position.x(), position.y(), position.z(), 0};
    G4double emComponents[6];

    G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldManager->GetDetectorField()->GetFieldValue(spacetimePoint, emComponents);
    G4double B[3] = {emComponents[0], emComponents[1], emComponents[2]};
    G4double normB = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

    G4double normMomentum = std::sqrt(pow(momentumDirection.x(), 2) + pow(momentumDirection.y(), 2) + pow(momentumDirection.z(), 2));
    G4double dotProd = (momentumDirection.x() * B[0]) + (momentumDirection.y() * B[1]) + (momentumDirection.z() * B[2]);

    G4double pitchAngleDeg = std::acos(dotProd / (normMomentum * normB)) * 180/3.14159265358979;

    // Write particle parameters to memory
    fRunAction->fBackscatteredParticleNames.push_back(particleName);
    fRunAction->fBackscatteredTrackWeights.push_back(trackWeight);
    fRunAction->fBackscatteredEnergieskeV.push_back(postStepKineticEnergy/keV);
    fRunAction->fBackscatteredPitchAnglesDeg.push_back(pitchAngleDeg);
    fRunAction->fBackscatterDirections.push_back({momentumDirection.x(), momentumDirection.y(), momentumDirection.z()});
    fRunAction->fBackscatterPositions.push_back({position.x()/m, position.y()/m, (position.z()/m) + 500000.0}); // Shift z-axis so we are writing altitude above sea level to file rather than the world coordinates

    // Kill particle
    track->SetTrackStatus(fStopAndKill);
  }



}

void SteppingAction::LogEnergy(G4int histogramAddress, G4double energy)
{
  // This is in a different function so the threadlock isn't in scope for all of every step.
  // Lock unlocks when it goes out of scope
  G4AutoLock lock(&aMutex); // Might not be necessary with thread-specific files
  fRunAction->fEnergyDepositionHistogram->AddCountToBin(histogramAddress, energy);
}