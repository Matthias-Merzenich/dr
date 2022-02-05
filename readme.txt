This repository contains Dean Hickerson's drifter search program for Conway's
Game of Life and life-like cellular automata. The initial commit contains the
program, documentation, and known rotors file exactly as I received them. Since
then, I have added this readme and made the following additional changes:

In dr.c:
  * Line 21: increased value of MAXKNOWN
  * Line 22: increased value of MAXFILESIZE
  * Line 83: changed file name to read knownrotors file from the current folder

In knownrotors:
  * Significantly expanded the list of known rotors


-Matthias Merzenich