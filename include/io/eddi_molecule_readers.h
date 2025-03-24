/**
 * @file eddi_molecule_readers.h
 * @author Luca Guffanti
 * @brief Contains function prototypes for reading molecules from different file formats 
 */

#ifndef __EDDI_MOLECULE_READERS_H__
#define __EDDI_MOLECULE_READERS_H__

#include "eddi_molecule.h"

/**
 * @brief Reads a PDB file containing molecular information and stores the molecule
 * in a molecule descriptor
 * 
 * @param filename input file name
 * @param molecule molecule object. Must be already initialized
 * @return true if the file is read completely
 * @return false if any problem occurs during the reading operation
 */
bool eddi_read_pdb(const char* filename, eddi_molecule_t* molecule);


#endif // __EDDI_MOLECULE_READERS_H__