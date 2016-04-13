/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Pennsylvania nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <sbpl/planners/mha_planner.h>
#include <sbpl/discrete_space_information/environment_roboarm_mha.h>

#include <iostream>
#include <memory>
#include <string>

#include <cmath>
#include <cstring>

#include <gflags/gflags.h>

DEFINE_string(stats_file, "", "File to dump stats in");

using namespace std;


/*******************************************************************************
 * planrobarm - An example of planning a robot arm with six degrees-of-freedom.
 *
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/robarm for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int planrobarm(char *envCfgFilename, char *islandsFilename,
               char *solFilename) {
  int b_ret = 0;
  double allocated_time_secs = 5.0; //in seconds
  MDPConfig MDPCfg;
  bool bforwardsearch = true;

  //Initialize Environment (should be called before initializing anything else)
  EnvironmentROBARMMHA environment_robarm;

  if (!environment_robarm.InitializeEnv(envCfgFilename, islandsFilename)) {
    printf("ERROR: InitializeEnv failed\n");
    throw new SBPL_Exception();
  }

  //Initialize MDP Info
  if (!environment_robarm.InitializeMDPCfg(&MDPCfg)) {
    printf("ERROR: InitializeMDPCfg failed\n");
    throw new SBPL_Exception();
  }

  // SBPLPlanner* planner = NULL;
  // planner = new ARAPlanner(&environment_robarm, bforwardsearch);

  MHAReplanParams planner_params(0.0);
  planner_params.use_lazy = false;
  // planner_params.meta_search_type = mha_planner::MetaSearchType::ROUND_ROBIN;
  planner_params.meta_search_type = mha_planner::MetaSearchType::DTS;
  planner_params.planner_type = mha_planner::PlannerType::SMHA;
  planner_params.mha_type = mha_planner::MHAType::PLUS;
  planner_params.inflation_eps = 5.0;
  planner_params.final_eps = planner_params.inflation_eps;
  planner_params.max_time = 10.0;
  planner_params.return_first_solution = true;
  planner_params.repair_time = -1;

  std::unique_ptr<MHAPlanner> planner;
  planner.reset(new MHAPlanner(&environment_robarm,
                               environment_robarm.GetNumHeuristics(), bforwardsearch));

  if (planner->set_start(MDPCfg.startstateid) == 0) {
    printf("ERROR: failed to set start state\n");
    throw new SBPL_Exception();
  }

  if (planner->set_goal(MDPCfg.goalstateid) == 0) {
    printf("ERROR: failed to set goal state\n");
    throw new SBPL_Exception();
  }

  printf("start planning...\n");
  vector<int> solution_stateIDs_V;
  int sol_cost;
  b_ret = planner->replan(&solution_stateIDs_V, planner_params, &sol_cost);
  // b_ret = planner->replan(allocated_time_secs, &solution_stateIDs_V);
  printf("done planning\n");
  std::cout << "size of solution=" << solution_stateIDs_V.size() << std::endl;

  FILE *fSol = fopen(solFilename, "w");

  if (fSol == NULL) {
    printf("ERROR: could not open solution file\n");
    throw new SBPL_Exception();
  }

  for (unsigned int i = 0; i < solution_stateIDs_V.size(); i++) {
    environment_robarm.PrintState(solution_stateIDs_V[i], true, fSol);
  }
  fclose(fSol);

  //print a path
  if (b_ret) {
    //print the solution
    printf("Solution is found\n");
  } else {
    printf("Solution does not exist\n");
  }

  // Write stats to file
  if (!FLAGS_stats_file.empty()) {
    FILE *fStats = fopen(FLAGS_stats_file.c_str(), "w");

    if (fStats == NULL) {
      printf("ERROR: could not open stats file\n");
      throw new SBPL_Exception();
    }

    if (b_ret) {
      vector<PlannerStats> planner_stats;
      planner->get_search_stats(&planner_stats);
      fprintf(fStats, "%f %d %d\n", planner_stats[0].time, planner_stats[0].expands,
              planner_stats[0].cost);
    } else {
      fprintf(fStats, "-1 -1 -1\n");
    }

    fclose(fStats);
  }

  fflush(NULL);

  return b_ret;
}

/*******************************************************************************
 * main - Parse command line arguments and launch one of the sbpl examples above.
 *        Launch sbpl with just the -h option for usage help.
 *
 * @param argc The number of command-line arguments
 * @param argv The command-line arguments
 *******************************************************************************/
int main(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Usage: ./test_roboarm <env_file> <islands_file> <solution_file>\n");
    exit(1);
  }

  google::ParseCommandLineFlags(&argc, &argv, true);
  planrobarm(argv[1], argv[2], argv[3]);
  return 0;
}
