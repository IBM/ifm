/*
 * FileIO/mapfile.C
 *
 * Authors: Lucas C. Villa Real (lucasvr@br.ibm.com)
 *
 * Description:
 *   Simple int:float map file parser
 *
 * (C) Copyright IBM Corp 2014
 */
#include "FileIO/mapfile.h"
#include <assert.h>

map<uint64_t,real_t> *MapFile::getMap(uint64_t fieldNumber, uint64_t *maxKey)
{
  assert(maxKey);

  ifstream infile(m_filename.c_str());
  string line;

  map<uint64_t,real_t> *m = new map<uint64_t,real_t>;
  real_t value = 0.0;
  uint64_t key = 0;

  while (getline(infile, line)) {
    if (line.size() == 0 || line[0] == '#')
      continue;
    if (line.find(":") != string::npos) {
      fprintf(stderr, "Syntax error on %s. A SPACE separated list of tokens is expected, as in:\n",
        m_filename.c_str());
      fprintf(stderr, "<Id> <ManningsN> <Retention>\n\n");
      delete m;
      return NULL;
    }
    istringstream ss(line);
    ss >> key;
    if (key > *maxKey)
      *maxKey = key;
    for (uint64_t i=0; i<fieldNumber; ++i)
      ss >> value;
    m->insert(pair<uint64_t,real_t>(key, value));
  }

  return m;
}

real_t *MapFile::getList(uint64_t fieldNumber, uint64_t *numValues)
{
  uint64_t maxKey = 0;
  map<uint64_t,real_t> *m = this->getMap(fieldNumber, &maxKey);
  if (m == NULL)
    return NULL;

  real_t *values = new real_t[maxKey+1];
  for (uint64_t i=0; i<maxKey+1; ++i)
    values[i] = 0.0;

  for (map<uint64_t,real_t>::iterator it=m->begin(); it != m->end(); ++it) {
    uint64_t key = (*it).first;
    real_t value = (*it).second;
    values[key] = value;
  }
  delete m;

  if (numValues)
    *numValues = maxKey+1;
  return values;
}
