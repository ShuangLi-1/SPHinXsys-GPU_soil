/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	base_body_part.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_BODY_PART_H
#define BASE_BODY_PART_H

#include "base_body.h"

#include <string>

namespace SPH
{
/**
 * @class BodyPart
 * @brief An auxiliary class for SPHBody to indicate a part of the body.
 */
using namespace std::placeholders;
class BodyPart
{
  protected:
    UniquePtrsKeeper<Entity> unique_variable_ptrs_;

  public:
    BodyPart(SPHBody &sph_body, const std::string &body_part_name);
    virtual ~BodyPart() {};
    SPHBody &getSPHBody() { return sph_body_; };
    SPHSystem &getSPHSystem() { return sph_body_.getSPHSystem(); };
    std::string getName() { return body_part_name_; };
    int getPartID() { return part_id_; };
    DiscreteVariable<UnsignedInt> *dvIndexList() { return dv_index_list_; };
    SingularVariable<UnsignedInt> *svRangeSize() { return sv_range_size_; };

  protected:
    SPHBody &sph_body_;
    int part_id_;
    std::string body_part_name_;
    BaseParticles &base_particles_;
    DiscreteVariable<UnsignedInt> *dv_index_list_;
    SingularVariable<UnsignedInt> *sv_range_size_;
    DiscreteVariable<int> *dv_body_part_indicator_;
    Vecd *pos_;
};

/**
 * @class BodyPartByParticle
 * @brief A body part with a collection of particles.
 */
class BodyPartByParticle : public BodyPart
{
  public:
    typedef BodyPartByParticle BaseIdentifier;
    IndexVector body_part_particles_; /**< Collection particle in this body part. */
    BaseParticles &getBaseParticles() { return base_particles_; };
    IndexVector &LoopRange() { return body_part_particles_; };
    size_t SizeOfLoopRange() { return body_part_particles_.size(); };
    BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name);
    virtual ~BodyPartByParticle() {};
    void setBodyPartBounds(BoundingBox bbox);
    BoundingBox getBodyPartBounds();

    template <typename SearchMethod>
    class TargetParticleMask : public SearchMethod
    {
      public:
        template <class ExecutionPolicy, typename EnclosureType, typename... Args>
        TargetParticleMask(ExecutionPolicy &ex_policy, EnclosureType &encloser, Args... args)
            : SearchMethod(std::forward<Args>(args)...), part_id_(encloser.part_id_),
              body_part_indicator_(encloser.dv_body_part_indicator_->DelegatedData(ex_policy)) {}
        virtual ~TargetParticleMask() {}

        bool isInRange(UnsignedInt index_i, UnsignedInt index_j)
        {
            return (body_part_indicator_[index_j] == part_id_) && SearchMethod::isInRange(index_i, index_j);
        }

      protected:
        int part_id_;
        int *body_part_indicator_;
    };
    
  protected:
    BoundingBox body_part_bounds_;
    bool body_part_bounds_set_;

    typedef std::function<bool(size_t)> TaggingParticleMethod;
    void tagParticles(TaggingParticleMethod &tagging_particle_method);
};

/**
 * @class BodyPartByCell
 * @brief A body part with a collection of cell lists.
 */
class BodyPartByCell : public BodyPart
{
  public:
    typedef BodyPartByCell BaseIdentifier;
    ConcurrentCellLists body_part_cells_; /**< Collection of cells to indicate the body part. */
    ConcurrentCellLists &LoopRange() { return body_part_cells_; };
    size_t SizeOfLoopRange();

    BodyPartByCell(RealBody &real_body, const std::string &body_part_name);
    virtual ~BodyPartByCell() {};
    DiscreteVariable<UnsignedInt> *getParticleIndex() { return dv_particle_index_; };
    DiscreteVariable<UnsignedInt> *getCellOffset() { return dv_cell_offset_; };

  protected:
    BaseCellLinkedList &cell_linked_list_;
    DiscreteVariable<UnsignedInt> *dv_particle_index_;
    DiscreteVariable<UnsignedInt> *dv_cell_offset_;
    typedef std::function<bool(Vecd, Real)> TaggingCellMethod;
    void tagCells(TaggingCellMethod &tagging_cell_method);
};

/**
 * @class BodyRegionByParticle
 * @brief A  body part with the collection of particles within by a prescribed shape.
 */
class BodyRegionByParticle : public BodyPartByParticle
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    BodyRegionByParticle(SPHBody &sph_body, Shape &body_part_shape);
    BodyRegionByParticle(SPHBody &sph_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByParticle() {};
    Shape &getBodyPartShape() { return body_part_shape_; };

  protected:
    Shape &body_part_shape_;
    bool tagByContain(size_t particle_index);
};

/**
 * @class BodySurface
 * @brief A  body part with the collection of particles at surface of a body
 */
class BodySurface : public BodyPartByParticle
{
  public:
    explicit BodySurface(SPHBody &sph_body);
    virtual ~BodySurface() {};

  protected:
    Real particle_spacing_min_;
    bool tagNearSurface(size_t particle_index);
};

/**
 * @class BodySurfaceLayer
 * @brief A  body part with the collection of particles within the surface layers of a body.
 */
class BodySurfaceLayer : public BodyPartByParticle
{
  public:
    explicit BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness = 3.0);
    virtual ~BodySurfaceLayer() {};

  private:
    Real thickness_threshold_;
    bool tagSurfaceLayer(size_t particle_index);
};

/**
 * @class BodyRegionByCell
 * @brief A body part with the cell lists within a prescribed shape.
 */
class BodyRegionByCell : public BodyPartByCell
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    BodyRegionByCell(RealBody &real_body, Shape &body_part_shape);
    BodyRegionByCell(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByCell() {};
    Shape &getBodyPartShape() { return body_part_shape_; };

  private:
    Shape &body_part_shape_;
    bool checkNotFar(Vecd cell_position, Real threshold);
};

/**
 * @class NearShapeSurface
 * @brief A body part with the cell lists near the surface of a prescribed shape.
 * @ details The body part shape can be that of the body,
 * or a sub shape of the body shape, or a shape independent of the body shape.
 * Note that only cells near the surface of the body part shape are included.
 */
class NearShapeSurface : public BodyPartByCell
{
  private:
    UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;

  public:
    NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    NearShapeSurface(RealBody &real_body, LevelSetShape &level_set_shape);
    explicit NearShapeSurface(RealBody &real_body);
    NearShapeSurface(RealBody &real_body, const std::string &sub_shape_name);
    virtual ~NearShapeSurface() {};
    LevelSetShape &getLevelSetShape() { return level_set_shape_; };

  private:
    LevelSetShape &level_set_shape_;
    bool checkNearSurface(Vecd cell_position, Real threshold);
};

class AlignedBoxPart
{
    UniquePtrKeeper<SingularVariable<AlignedBox>> sv_aligned_box_keeper_;

  public:
    AlignedBoxPart(const std::string &name, const AlignedBox &aligned_box);
    virtual ~AlignedBoxPart() {};
    SingularVariable<AlignedBox> *svAlignedBox() { return sv_aligned_box_keeper_.getPtr(); };
    AlignedBox &getAlignedBox() { return aligned_box_; };

  protected:
    AlignedBox &aligned_box_;
};

class AlignedBoxPartByParticle : public BodyPartByParticle, public AlignedBoxPart
{
  public:
    AlignedBoxPartByParticle(RealBody &real_body, const AlignedBox &aligned_box);
    virtual ~AlignedBoxPartByParticle() {};

  protected:
    bool tagByContain(size_t particle_index);
};

class AlignedBoxPartByCell : public BodyPartByCell, public AlignedBoxPart
{
  public:
    AlignedBoxPartByCell(RealBody &real_body, const AlignedBox &aligned_box);
    virtual ~AlignedBoxPartByCell() {};

  protected:
    bool checkNotFar(Vecd cell_position, Real threshold);
};
} // namespace SPH
#endif // BASE_BODY_PART_H
