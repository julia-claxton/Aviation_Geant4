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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Creates particles that will be propagated in the simulation

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include <math.h>


PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fPrimaryMessenger(0),
  fBeamEnergy(100.),
  fBeamPitchAngle(40.),
  fInitialParticleAlt(450.0),
  fPI(3.14159265359),
  fRad2Deg(180.0 / 3.14159265359),
  fSourceType("e-")
{
  fPrimaryMessenger = new PrimaryGeneratorMessenger(this);
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Select input particle type
  fParticleGun  = new G4ParticleGun();
  G4ParticleDefinition* inputParticle = G4ParticleTable::GetParticleTable()->FindParticle(fSourceType); // Electron = "e-", proton = "proton", photon = "gamma"
  fParticleGun->SetParticleDefinition(inputParticle);
  
  // Set up container for particle properties
  ParticleSample* r = new ParticleSample();

  // Set particle energy
  r->energy = fBeamEnergy * keV;
 
  // Set starting location
  r->xPos = 0; 
  r->yPos = 0;
  r->zPos = (fInitialParticleAlt - 500.0)*km; // Subtraction due to coordinate axis location in middle of world volume

  /*
  // Set velocity. Starts electrons with gyromotion about field line at a given pitch angle.
  G4double pitchAngle = fBeamPitchAngle * 3.14159265359 / 180.0; // Convert degrees to radians
  G4double gyroPhase  = G4UniformRand() * 2. * fPI;

  // Initial momentum direction of particles in local frame of magnetic field (i.e B is parallel to -z)
  double vx0 = std::sin(pitchAngle)*std::cos(gyroPhase);
  double vy0 = std::sin(pitchAngle)*std::sin(gyroPhase);
  double vz0 = -std::cos(pitchAngle);
  
  // Now we need to rotate out of inclined B-field frame into the world frame.
  // B lies in the y-z plane, meaning we need to rotate about the x-axis to get into the world coordinates.
  
  // Get B vector
  G4double spacetimePoint[4] = {r->xPos, r->xPos, r->zPos, 0};
  G4double emComponents[6];

  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldManager->GetDetectorField()->GetFieldValue(spacetimePoint, emComponents);
  G4double B[3] = {emComponents[0], emComponents[1], emComponents[2]};
  G4double normB = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  
  // Get necessary tilt angle
  G4double tilt_angle_rad = std::acos(std::abs(B[2]/normB));

  // Rotate coordinates into world frame
  r->xDir = vx0;
  r->yDir = (std::cos(tilt_angle_rad) * vy0) - (std::sin(tilt_angle_rad) * vz0);
  r->zDir = (std::sin(tilt_angle_rad) * vy0) + (std::cos(tilt_angle_rad) * vz0);

  // Verify that pitch angle generation is correct
  double normMomentum = std::sqrt(pow(r->xDir, 2) + pow(r->yDir, 2) + pow(r->zDir, 2));
  double dotProd = (r->xDir * B[0]) + (r->yDir * B[1]) + (r->zDir * B[2]);
  double generatedPitchAngle_deg = std::acos(dotProd / (normMomentum * normB)) * 180/3.14159265358979;

  if( abs(generatedPitchAngle_deg - fBeamPitchAngle) > 1){
    G4cout << "** ERROR: Primary generated with pitch angle >1º different than user-specified pitch angle. You should never see this. Please email julia.claxton@colorado.edu with this error and the conditions that produced it." << G4endl;
    throw;
  }
  */

  // Send particles straight down: perpendicular incidence
  r->xDir = 0;
  r->yDir = 0;
  r->zDir = -1;

  // Communicate to particle gun
  fParticleGun->SetParticlePosition(G4ThreeVector(r->xPos, r->yPos, r->zPos)); 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(r->xDir, r->yDir, r->zDir));
  fParticleGun->SetParticleEnergy(r->energy);
  
  // Geant method to create initial particle with the above properties 
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Free memory from ParticleSample struct
  delete r;
}


