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

#ifndef __ENVIRONMENT_MHA_H_
#define __ENVIRONMENT_MHA_H_

#include <sbpl/discrete_space_information/environment.h>
#include <sbpl/utils/utils.h>

class EnvironmentMHA : public DiscreteSpaceInformation
{
  public:
    /**
     * \brief heuristic estimate from state FromStateID to state ToStateID
     */
    virtual int GetFromToHeuristic(int q_id, int FromStateID, int ToStateID)
    {
      return 0;
    }

    /**
     * \brief heuristic estimate from state with stateID to goal state
     */
    virtual int GetGoalHeuristic(int q_id, int stateID) = 0;

    /**
     * \brief heuristic estimate from start state to state with stateID
     */
    virtual int GetStartHeuristic(int q_id, int stateID) 
    {
      return 0;
    }

    /**
     * \brief GetSuccs methods that inform the environment which queue the states are being expanded from
     */
    using DiscreteSpaceInformation::GetSuccs;
    using DiscreteSpaceInformation::GetLazySuccs;
    using DiscreteSpaceInformation::GetLazyPreds;
    using DiscreteSpaceInformation::GetTrueCost;

    virtual void GetSuccs(int q_id, int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV) 
    {
      SBPL_ERROR("ERROR: GetSuccs with q_id is not implemented for this environment!\n");
      throw new SBPL_Exception();
    }

    virtual void GetLazySuccs(int q_id, int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV, std::vector<bool>* isTrueCost){
      SBPL_ERROR("ERROR: GetLazySuccs with q_id is not implemented for this environment!\n");
      throw new SBPL_Exception();
    };

    virtual void GetLazyPreds(int q_id, int SourceStateID, std::vector<int>* SuccIDV, std::vector<int>* CostV, std::vector<bool>* isTrueCost){
      SBPL_ERROR("ERROR: GetLazyPreds with q_id is not implemented for this environment!\n");
      throw new SBPL_Exception();
    };

    virtual int GetTrueCost(int q_id, int parentID, int childID){
      SBPL_ERROR("ERROR: GetTrueCost with q_id is not implemented for this environment!\n");
      throw new SBPL_Exception();
      return -1;
    };
};

#endif

