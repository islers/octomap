/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
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

#ifndef OCTOMAP_OCCUPANCY_OCTREE_STEREO_H
#define OCTOMAP_OCCUPANCY_OCTREE_STEREO_H


#include <vector>
#include "OccupancyOcTreeBase.h"


namespace octomap {

  /**
   * Stereo sensor model implementation for Occupancy Octrees (e.g. for mapping).
   * Uses a quadratic inverse depth model for computing logodds updates.
   * AbstractOccupancyOcTree serves as a common
   * base interface for all these classes.
   * Each class used as NODE type needs to be derived from
   * OccupancyOcTreeNode.
   *
   * This tree implementation has a maximum depth of 16. 
   * At a resolution of 1 cm, values have to be < +/- 327.68 meters (2^15)
   *
   * This limitation enables the use of an efficient key generation 
   * method which uses the binary representation of the data.
   *
   * \note The tree does not save individual points.
   *
   * \tparam NODE Node class to be used in tree (usually derived from
   *    OcTreeDataNode)
   */
  template <class NODE>
  class OccupancyOcTreeStereo : public OccupancyOcTreeBase<NODE> {

  public:
    /// Default constructor, sets resolution of leafs
    OccupancyOcTreeStereo(double resolution, double coeff, double max_range);
    virtual ~OccupancyOcTreeStereo();

    /// Copy constructor
    OccupancyOcTreeStereo(const OccupancyOcTreeStereo<NODE>& rhs);

    /**
    * Integrate a Pointcloud (in global reference frame), parallelized with OpenMP.
    * Special care is taken that each voxel
    * in the map is updated only once, and occupied nodes have a preference over free ones.
    * This avoids holes in the floor from mutual deletion and is more efficient than the plain
    * ray insertion in insertPointCloudRays().
    *
    * @note replaces insertScan()
    *
    * @param scan Pointcloud (measurement endpoints), in global reference frame
    * @param sensor_origin measurement origin in global reference frame
    * @param sensor_orientation measurement orientation as vector pointing in the direction of 
    *   the camera center in global reference frame
    * @param maxrange maximum range for how long individual beams are inserted 
    *   (default -1: use max_range set with constructor)
    * @param lazy_eval whether update of inner nodes is omitted after the update (default: false).
    *   This speeds up the insertion, but you need to call updateInnerOccupancy() when done.
    * @param discretize whether the scan is discretized first into octree key cells (default: false).
    *   This reduces the number of raycasts using computeDiscreteUpdate(), resulting in a potential speedup.*
    */
    virtual void insertPointCloud(const Pointcloud& scan, const octomap::point3d& sensor_origin,
                                  const octomath::Vector3& sensor_orientation, double maxrange=-1., 
                                  bool lazy_eval = false, bool discretize = false);

    /**
     * Integrate occupancy measurement with dependence on distance from sensor to measurement.
     *
     * @param key OcTreeKey of the NODE that is to be updated
     * @param occupied true if the node was measured occupied, else false
     * @param origin = sensor origin in world frame
     * @return pointer to the updated NODE
     */
    virtual NODE* updateNode(const OcTreeKey& key, bool occupied, float d);

    /** returns the distance dependent logodds update value for a measurement at a certain depth.
     * @param occupied true if the node was measured occupied, else false
     * @param distance distance of measurement from sensor origin
     */
    float getLogOdds( bool occupied, float distance );

    inline double getCoeff() const { return coefficient; }
    inline double getMaxRange() const { return maximum_range; }
    inline void setCoeff(double coeff) { coefficient = coeff; }
    inline void setMaxRange(double max_range) { maximum_range = max_range; }

  protected:
    void generateLUTs();

    double coefficient; // Stereo coefficient
    double maximum_range; // Maximum sensor range
    std::vector<float> hit_LUT;
    std::vector<float> miss_LUT;
  };

} // namespace

#include "octomap/OccupancyOcTreeStereo.hxx"

#endif
