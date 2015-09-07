/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Stereod on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
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

namespace octomap {

  template <class NODE>
  OccupancyOcTreeStereo<NODE>::OccupancyOcTreeStereo(double resolution, double coeff, double max_range)
    : OccupancyOcTreeBase<NODE>(resolution), coefficient(coeff), maximum_range(max_range) 
  {
    generateLUTs();
  }
  
  template <class NODE>
  OccupancyOcTreeStereo<NODE>::~OccupancyOcTreeStereo(){
  }

  template <class NODE>
  OccupancyOcTreeStereo<NODE>::OccupancyOcTreeStereo(const OccupancyOcTreeStereo<NODE>& rhs) :
      OcTreeBaseImpl<NODE,AbstractOccupancyOcTree>(rhs),  
      coefficient(rhs.coefficient), maximum_range(rhs.maximum_range),
      hit_LUT(rhs.hit_LUT), miss_LUT(rhs.miss_LUT)
  {
    this->use_bbx_limit = rhs.use_bbx_limit,
    this->bbx_min = rhs.bbx_min; 
    this->bbx_max = rhs.bbx_max;
    this->bbx_min_key = rhs.bbx_min_key;
    this->bbx_max_key = rhs.bbx_max_key;
    this->use_change_detection = rhs.use_change_detection; 
    this->changed_keys = rhs.changed_keys;
    this->clamping_thres_min = rhs.clamping_thres_min;
    this->clamping_thres_max = rhs.clamping_thres_max;
    this->prob_hit_log = rhs.prob_hit_log;
    this->prob_miss_log = rhs.prob_miss_log;
    this->occ_prob_thres_log = rhs.occ_prob_thres_log;
  }

  template <class NODE>
  void OccupancyOcTreeStereo<NODE>::insertPointCloud(const Pointcloud& scan, 
                                                     const octomap::point3d& sensor_origin,
                                                     const octomath::Vector3& sensor_orientation,
                                                     double maxrange, bool lazy_eval, bool discretize) {
    if(maxrange == -1.0)
      maxrange = maximum_range;
    else if(maxrange > maximum_range)
      maxrange = maximum_range;

    octomath::Vector3 orient = sensor_orientation.normalized();

    KeySet free_cells, occupied_cells;
    
    if (discretize)
      OccupancyOcTreeBase<NODE>::computeDiscreteUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);
    else
      OccupancyOcTreeBase<NODE>::computeUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);

    // insert data into tree  -----------------------
    for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
      float d = fabs(((OccupancyOcTreeBase<NODE>::keyToCoord(*it)) - sensor_origin).dot(orient));
      updateNode(*it, false, d);
    }
    for (KeySet::iterator it = occupied_cells.begin(); it != occupied_cells.end(); ++it) {
      float d = fabs(((OccupancyOcTreeBase<NODE>::keyToCoord(*it)) - sensor_origin).dot(orient));
      updateNode(*it, true, d);
    }
  }

  template <class NODE>
  NODE* OccupancyOcTreeStereo<NODE>::updateNode(const OcTreeKey& key, 
                                                bool occupied, 
                                                float d) {
     float logOdds = getLogOdds(occupied,d);
     return OccupancyOcTreeBase<NODE>::updateNode(key, logOdds, false);
  }
  
  template <class NODE>
  float OccupancyOcTreeStereo<NODE>::getLogOdds( bool occupied, float d )
  {
    float logOdds;
     if (d > maximum_range) d = maximum_range;
     if (d < 0.0) d = 0.0;
     int lutIndex = floor(d/this->resolution);
     if (occupied)
       logOdds = hit_LUT.at(lutIndex);
     else
       logOdds = miss_LUT.at(lutIndex);
     
     return logOdds;
  }

  template <class NODE>
  void OccupancyOcTreeStereo<NODE>::generateLUTs()
  {
    int levels = floor(maximum_range/this->resolution);
    for(int i = 0; i <= levels; i++)
    {
      float factor = erf(coefficient/(i*i*this->resolution*this->resolution));
      float hitLogOdds = logodds(0.5 + ((probability(this->prob_hit_log) - 0.5)*factor));
      float missLogOdds = logodds(0.5 + ((probability(this->prob_miss_log) - 0.5)*factor));
      hit_LUT.push_back(hitLogOdds);
      miss_LUT.push_back(missLogOdds);
    }
  }

} // namespace
