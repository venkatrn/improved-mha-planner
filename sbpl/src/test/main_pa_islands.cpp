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
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// #include <sbpl/headers.h>
#include <sbpl/planners/mha_planner.h>
#include <sbpl/discrete_space_information/environment_navxythetalat_mha.h>

#include <gflags/gflags.h>

DEFINE_bool(visualization, true, "Visualize states expanded during search (slows down the search)");

using namespace std;

/*******************************************************************************
 * planxythetalat
 * @brief An example of planning in three-dimensional space (x,y,theta)
 *
 * @param plannerType The type of planner to be used in this example
 * @param envCfgFilename The environment config file. See
 *                       sbpl/env_examples/nav3d/ for examples
 * @param motPrimFilename The motion primitives file. See
 *                        sbpl/matlab/mprim/ for examples
 * @return 1 if the planner successfully found a solution; 0 otherwise
 *******************************************************************************/
int planxythetalat(char* envCfgFilename, char* motPrimFilename, char* islandsFilename, char* solFilename)
{
    MDPConfig MDPCfg;

    MHAReplanParams planner_params(0.0);
    planner_params.use_lazy = false;
    // planner_params.meta_search_type = mha_planner::MetaSearchType::ROUND_ROBIN;
    planner_params.meta_search_type = mha_planner::MetaSearchType::DTS;
    planner_params.planner_type = mha_planner::PlannerType::SMHA;
    planner_params.mha_type = mha_planner::MHAType::PLUS;
    planner_params.final_eps = planner_params.inflation_eps;
    planner_params.inflation_eps = 500.0;
    planner_params.max_time = 10.0;
    planner_params.return_first_solution = true;
    planner_params.repair_time = -1;

    // set the perimeter of the robot (it is given with 0,0,0 robot ref. point for which planning is done)
    vector<sbpl_2Dpt_t> perimeterptsV;
    sbpl_2Dpt_t pt_m;
    double halfwidth = 0.01; //0.3;
    double halflength = 0.01; //0.45;
    pt_m.x = -halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = -halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);
    pt_m.x = -halflength;
    pt_m.y = halfwidth;
    perimeterptsV.push_back(pt_m);

    // clear the footprint
    perimeterptsV.clear();

    // Initialize Environment (should be called before initializing anything else)
    EnvironmentNAVXYTHETALAT environment_navxythetalat;
		environment_navxythetalat.SetVisualization(FLAGS_visualization);

    if (!environment_navxythetalat.InitializeEnv(envCfgFilename, perimeterptsV, motPrimFilename, islandsFilename)) {
        printf("ERROR: InitializeEnv failed\n");
        throw new SBPL_Exception();
    }

    // Initialize MDP Info
    if (!environment_navxythetalat.InitializeMDPCfg(&MDPCfg)) {
        printf("ERROR: InitializeMDPCfg failed\n");
        throw new SBPL_Exception();
    }

    // plan a path

    unique_ptr<MHAPlanner> planner;
    int num_heuristics = environment_navxythetalat.GetNumHeuristics();
    planner.reset(new MHAPlanner(&environment_navxythetalat, num_heuristics, true));


    // set planner properties
    if (planner->set_start(MDPCfg.startstateid) == 0) {
        printf("ERROR: failed to set start state\n");
        throw new SBPL_Exception();
    }
    if (planner->set_goal(MDPCfg.goalstateid) == 0) {
        printf("ERROR: failed to set goal state\n");
        throw new SBPL_Exception();
    }
    // plan
    vector<int> solution_stateIDs_V;
    int sol_cost;
    bool b_ret;
    printf("start planning...\n");
    b_ret = planner->replan(&solution_stateIDs_V, planner_params, &sol_cost);
    printf("done planning\n");
    printf("size of solution=%d\n", (unsigned int)solution_stateIDs_V.size());

    environment_navxythetalat.PrintTimeStat(stdout);

    // write solution to sol.txt
    FILE* fSol = fopen(solFilename, "w");
    if (fSol == NULL) {
        printf("ERROR: could not open solution file\n");
        throw new SBPL_Exception();
    }

    // write the discrete solution to file
    //	for (size_t i = 0; i < solution_stateIDs_V.size(); i++) {
    //		int x;
    //		int y;
    //		int theta;
    //		environment_navxythetalat.GetCoordFromState(solution_stateIDs_V[i], x, y, theta);
    //
    //		fprintf(fSol, "%d %d %d\t\t%.3f %.3f %.3f\n", x, y, theta,
    //              DISCXY2CONT(x, 0.1), DISCXY2CONT(y, 0.1), DiscTheta2Cont(theta, 16));
    //	}

    // write the continuous solution to file
    vector<sbpl_xy_theta_pt_t> xythetaPath;
    environment_navxythetalat.ConvertStateIDPathintoXYThetaPath(&solution_stateIDs_V, &xythetaPath);
    printf("solution size=%d\n", (unsigned int)xythetaPath.size());
    for (unsigned int i = 0; i < xythetaPath.size(); i++) {
        fprintf(fSol, "%.3f %.3f %.3f\n", xythetaPath.at(i).x, xythetaPath.at(i).y, xythetaPath.at(i).theta);
    }
    fclose(fSol);

    environment_navxythetalat.PrintTimeStat(stdout);

    // print a path
    if (b_ret) {
        // print the solution
        printf("Solution is found\n");
    }
    else {
        printf("Solution does not exist\n");
    }

    fflush(NULL);
    return b_ret;
}

int main(int argc, char *argv[])
{
    if (argc < 5) {
      printf("Usage: ./main_pa_islands <env_file> <mprim_file> <islands_file> <solution_file>\n");
      exit(1);
    }
		google::ParseCommandLineFlags(&argc, &argv, true);
    char *env_file = argv[1];
    char *mprim_file = argv[2];
    char *islands_file = argv[3];
    char *sol_file = argv[4];
    planxythetalat(env_file, mprim_file, islands_file, sol_file);
    return 0;
}
