#include <sbpl/discrete_space_information/environment_roboarm_mha.h>

#include <string>
#include <vector>

#include <cstring>

using namespace std;

bool EnvironmentROBARMMHA::InitializeEnv(const char *sEnvFile, const char *sIslandsFile) {
  InitializeEnv(sEnvFile);
  FILE *fIslands = fopen(sIslandsFile, "r");

  if (fIslands == NULL) {
    SBPL_ERROR("ERROR: unable to open %s\n", sIslandsFile);
    throw new SBPL_Exception();
  }

  const bool success = ReadIslands(fIslands);
  fclose(fIslands);
  return success;
}

int EnvironmentROBARMMHA::GetGoalHeuristic(int q_id, int stateID) {
#if USE_HEUR==0
    return 0;
#endif
  if (q_id == 0) {
  return GetGoalHeuristic(stateID);
  }

#if DEBUG
    if (stateID >= (int)EnvROBARM.StateID2CoordTable.size()) {
        SBPL_ERROR("ERROR in EnvROBARM... function: stateID illegal\n");
        throw new SBPL_Exception();
    }
#endif

    if (q_id > static_cast<int>(islands_.size())) {
      SBPL_ERROR("Requested q_id is greater than number of islands\n");
    }

    //define this function if it used in the planner (heuristic forward search would use it)
    vector<double> island = islands_[q_id - 1];

    double angles[NUMOFLINKS];
    short unsigned int coord[NUMOFLINKS];
    short unsigned int endeffx, endeffy;
    ComputeCoord(island.data(), coord);
    ComputeContAngles(coord, angles);
    ComputeEndEffectorPos(angles, &endeffx, &endeffy);

    //create the start state
    EnvROBARMHashEntry_t* HashEntry = CreateNewHashEntry(coord, NUMOFLINKS, endeffx, endeffy);
    return EnvironmentROBARM::GetFromToHeuristic(stateID, HashEntry->stateID);
}

void EnvironmentROBARMMHA::GetSuccs(int q_id, int SourceStateID,
                                    std::vector<int> *SuccIDV, std::vector<int> *CostV) {
  GetSuccs(SourceStateID, SuccIDV, CostV);
}

bool EnvironmentROBARMMHA::ReadIslands(FILE *fIslands) {
  char sTemp[1024], sExpected[1024];
  int numIslands = 0;

  SBPL_PRINTF("Reading in islands...");

  //read in the number of islands
  strcpy(sExpected, "num_islands:");

  if (fscanf(fIslands, "%s", sTemp) == 0) {
    return false;
  }

  if (strcmp(sTemp, sExpected) != 0) {
    SBPL_ERROR("ERROR: expected %s but got %s\n", sExpected, sTemp);
    return false;
  }

  if (fscanf(fIslands, "%d", &numIslands) == 0) {
    return false;
  }

  if (numIslands < 0) {
    SBPL_ERROR("ERROR: number of islands should be non-negative\n");
    return false;
  }

  islands_.reserve(numIslands);

  //read in the islands
  for (int ii = 0; ii < numIslands; ++ii) {
    vector<double> joint_angles(NUMOFLINKS, 0);

    for (int jj = 0; jj < NUMOFLINKS; ++jj) {
      if (fscanf(fIslands, "%lf", &joint_angles[jj]) != 1) {
        SBPL_ERROR("ERROR: incorrect format for islands file\n");
        return false;
      }
    }

    // TODO: validate islands against environment bounds and obstacles.
    islands_.push_back(joint_angles);
  }

  SBPL_PRINTF("done");
  SBPL_PRINTF("Islands:");
  for (int ii = 0; ii < numIslands; ++ii) {
    string str = "";

    for (int jj = 0; jj < NUMOFLINKS; ++jj) {
      str += to_string(islands_[ii][jj]) + " ";
    }

    SBPL_PRINTF("%s", str.c_str());
  }

  return true;
}
