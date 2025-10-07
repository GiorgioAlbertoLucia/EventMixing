#pragma once 

#include <TTree.h>

/**
 * Abstract version of the candidate to mix
*/
class Candidate 
{
    public:
        virtual ~Candidate() = default;
        virtual void setBranchAddress(TTree * tree) = 0;
};

/**
 * Structure to define the brackets of hadrons in a given collision (indices of hadrons produced in the same collision)
*/
struct CollHadBracket
{
    int CollID, fHadStartIndex = -1, fHadEndIndex = -1;
    void SetMin(int min) {
        fHadStartIndex = min;
    }
    void SetMax(int max) {
        fHadEndIndex = max;
    }
    int GetMin() {
        return fHadStartIndex;
    }
    int GetMax() {
        return fHadEndIndex;
    }
};
