

#include "SteppingActionMessenger.hh"

#include "SteppingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


SteppingActionMessenger::SteppingActionMessenger(SteppingAction* step)
  :G4UImessenger(),
   fSteppingAction(step) 
{
  fSteppingActionDir = new G4UIdirectory("/dataCollection/");
  fSteppingActionDir->SetGuidance("Select data collection type.");

  fcmd1 = new G4UIcmdWithAString("/dataCollection/setBackscatterFilename", this);
  fcmd1->SetParameterName("Enter file name.", true);
  fcmd1->SetDefaultValue("backscatter.csv");

  fcmd2 = new G4UIcmdWithADouble("/dataCollection/setCollectionAltitude",this);
  fcmd2->SetParameterName("Altitude to record upgoing particles at",true);
  fcmd2->SetDefaultValue(450.0);
}

SteppingActionMessenger::~SteppingActionMessenger()
{
  delete fSteppingActionDir;
  delete fcmd1;
  delete fcmd2;
}

void SteppingActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fcmd1){fSteppingAction->SetBackscatterFilename(newValue);}
  if(command == fcmd2){fSteppingAction->SetCollectionAltitude(std::stod(newValue));}
}
