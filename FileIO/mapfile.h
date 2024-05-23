/*
 * FileIO/mapfile.h
 *
 * Authors: Lucas C. Villa Real (lucasvr@br.ibm.com)
 *
 * Description:
 *   Simple int:float map file parser
 *
 * (C) Copyright IBM Corp 2014
 */
#ifndef __MAPFILE_H
#define __MAPFILE_H

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include "ifm_common.h"

using namespace std;

class MapFile {
  public:
    MapFile(string filename) : m_filename(filename) { };
    ~MapFile() { };

    // Tells how many entries (or mappings) we have on the file.
    //int  numEntries();

    // Tells how many fields a given rule has.
    //int  numFields();

    // Returns an <int,real_t> map where <real_t> is the rule's value at fieldNumber.
    // fieldNumber starts at index 1.
    // maxKey contains the maximum value assigned to <int>.
    map<uint64_t,real_t> *getMap(uint64_t fieldNumber, uint64_t *maxKey);

    // Returns a real_t array with the rule's values at fieldNumber.
    // fieldNumber starts at index 1.
    // numValues contains the number of entries in the returned array.
    real_t *getList(uint64_t fieldNumber, uint64_t *numValues=NULL);

  private:
    string m_filename;
};

#endif /* __MAPFILE_H */
