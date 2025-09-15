# Geant4_Atmospheric_EPP

Forked from Grant Berland. Readme TBD until code is finished.

need to edit `EDIT_THIS_FILE.mac` with your build directory. all user-editable simulation parameters also live there so take a look before running

How to build executable:
... TODO

```
cd "path/to/g4epp_source"

cmake -DCMAKE_INSTALL_PREFIX="/path/to/geant4-install" -DGeant4_DIR="/path/to/geant4-install/lib" "path/to/g4epp_source"

make
```


How to run executable:

`./G4EPP <number of particles> <particle> <energy> <pitch angle>`

particle:
* e- = electrons
* proton = protons
* gamma = photon
* For other particles, see: https://fismed.ciemat.es/GAMOS/GAMOS_doc/GAMOS.5.1.0/x11519.html

energy in kev

pitch angle in degrees

**section**
julia script is for batch simulation to create data for G4EPP_2.0 (link)

add columns of msis file