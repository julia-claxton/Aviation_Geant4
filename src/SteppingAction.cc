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
  const G4double preStepKineticEnergy = step->GetPreStepPoint()->GetKineticEnergy();
  const G4double postStepKineticEnergy = step->GetPostStepPoint()->GetKineticEnergy();
  const G4double trackWeight = track->GetWeight();

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
  // Energy Spectrum Tracking
  // ===========================
  G4double zPos = position.z(); // Particle altitude in world coordinates
  G4double alt_km = (zPos/km) + 500.0; // Altitude from ground level. + 500 because z = 0 in the simulation is at 500 km above ground level.

  G4double sampleAltitudes_km[] = {200.0, 110.0, 90.0, 400.0, 70.0}; // Altitudes to sample the energy spectrum at, km

  G4double preStepAlt_km  = (step->GetPreStepPoint()->GetPosition().z()/km) + 500.0;
  G4double postStepAlt_km = (step->GetPostStepPoint()->GetPosition().z()/km) + 500.0;

  // Loop over every sample altitude and see if this particle has crossed that plane
  for(int i = 0; i < std::end(sampleAltitudes_km) - std::begin(sampleAltitudes_km); i++){
    bool crossedPlane = (preStepAlt_km > sampleAltitudes_km[i]) != (postStepAlt_km > sampleAltitudes_km[i]);
    if( crossedPlane == false ){continue;} // Don't proceed if we haven't crossed the plane
  
    G4double crossingEnergy = postStepKineticEnergy;

    // If it's not safe to use pre- and post- kinetic energy interchangably, do an interpolation to get approximate energy at the plane crossing.
    if( std::abs((postStepKineticEnergy - preStepKineticEnergy)/preStepKineticEnergy) > 1e-10 ){
      G4double t = (sampleAltitudes_km[i] - preStepAlt_km) / (postStepAlt_km - preStepAlt_km);
      crossingEnergy = preStepKineticEnergy + (t * (postStepKineticEnergy - preStepKineticEnergy)); // Linear interpolation
      
      G4cout << "WARNING: Interpolated crossing energy." << G4endl; // Warn about it
    }

    G4cout << particleName << " crossed plane at " << sampleAltitudes_km[i] << " km with " << crossingEnergy/keV << " keV" << G4endl;
    // TODO record into histogram
    // TODO separate energy spectra for p+ e- and alpha

  }
}

void SteppingAction::LogEnergy(G4int histogramAddress, G4double energy)
{
  // This is in a different function so the threadlock isn't in scope for all of every step.
  // Lock unlocks when it goes out of scope
  G4AutoLock lock(&aMutex); // Might not be necessary with thread-specific files
  //fRunAction->fEnergyDepositionHistogram->AddCountToBin(histogramAddress, energy);
}