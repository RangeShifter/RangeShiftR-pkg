#ifndef TRAITTH
#define TRAITTH

#include <vector>
#include<set>

#include "parameters.h"
#include "Mutation.h"
#include "ProtoTrait.h"

using namespace std;

class TTrait {
public:
    /** Mutation procedure, perform mutations on the genes sequence. **/
    virtual void mutate() = 0;

    virtual unique_ptr<TTrait> clone() const = 0; //copies parameters (if not static) not gene seqeunces

    /** Inheritance procedure, creates a new trait from mother's and father's traits
  * @param mother the mother's trait
  * @param father the father's trait
  **/
    virtual void inherit(TTrait*, set<unsigned int> const& , sex_t, int ) = 0;

    /** sequence accessor.
  * @return the sequence pointer
  **/
  //  virtual   void** get_sequence() const = 0;
    virtual int getNLoci() const = 0;
    virtual float getMutationRate() const = 0;
    virtual bool isInherited() const = 0;
    //virtual   vector<vector<shared_ptr<Mutation>>> get_mutations() const = 0;
    //virtual   vector<vector<char>> get_mutations() const = 0;
    virtual float getSelectionCoefAtLoci(short chromosome, int i) const = 0;
    virtual int countHeterozygoteLoci() const = 0;
    virtual bool isHeterozygoteAtLoci(int loci) const = 0;
    virtual float express() = 0;
    virtual ~TTrait() { }
};

extern RSrandom* pRandom;

#endif