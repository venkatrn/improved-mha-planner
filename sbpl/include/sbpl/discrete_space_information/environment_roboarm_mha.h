#pragma once

#include <sbpl/discrete_space_information/environment_mha.h>
#include <sbpl/discrete_space_information/environment_robarm.h>

#include <vector>

class EnvironmentROBARMMHA: public EnvironmentMHA, public EnvironmentROBARM
{
  public:
    virtual bool InitializeEnv(const char* sEnvFile, const char* sIslandsFile);
    virtual int GetGoalHeuristic(int q_id, int stateID);
    virtual void GetSuccs(int q_id, int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV);
    int GetNumHeuristics() {
      // Num islands + 1 admissible heuristic
      return static_cast<int>(islands_.size()) + 1;
    }

    using EnvironmentROBARM::GetSuccs;
    using EnvironmentROBARM::GetGoalHeuristic;
    using EnvironmentROBARM::InitializeEnv;

  private:
    // Each island is a double vector of length NUMOFLINKS (continuous angles)
    std::vector<std::vector<double>> islands_;
    bool ReadIslands(FILE *fIslands);
};
