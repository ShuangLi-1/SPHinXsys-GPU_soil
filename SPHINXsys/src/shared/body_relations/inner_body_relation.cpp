/**
 * @file 	inner_body_relation.cpp
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "inner_body_relation.h"

#include "base_particle_dynamics.h"
#include "cell_linked_list.hpp"

namespace SPH
{
	//=================================================================================================//
	InnerRelation::InnerRelation(RealBody &real_body)
		: BaseInnerRelation(real_body), get_inner_neighbor_(real_body),
		  cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
	//=================================================================================================//
	void InnerRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		cell_linked_list_.searchNeighborsByParticles(
			sph_body_, inner_configuration_,
			get_single_search_depth_, get_inner_neighbor_);
	}
	//=================================================================================================//
	AdaptiveInnerRelation::
		AdaptiveInnerRelation(RealBody &real_body)
		: BaseInnerRelation(real_body), total_levels_(0),
		  get_adaptive_inner_neighbor_(real_body)
	{
		MultilevelCellLinkedList &multi_level_cell_linked_list =
			DynamicCast<MultilevelCellLinkedList>(this, real_body.getCellLinkedList());
		cell_linked_list_levels_ = multi_level_cell_linked_list.getMeshLevels();
		total_levels_ = cell_linked_list_levels_.size();
		for (size_t l = 0; l != total_levels_; ++l)
		{
			get_multi_level_search_depth_.push_back(
				adaptive_search_depth_ptr_vector_keeper_
					.createPtr<SearchDepthAdaptive>(real_body, cell_linked_list_levels_[l]));
		}
	}
	//=================================================================================================//
	void AdaptiveInnerRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		for (size_t l = 0; l != total_levels_; ++l)
		{
			cell_linked_list_levels_[l]->searchNeighborsByParticles(
				sph_body_, inner_configuration_,
				*get_multi_level_search_depth_[l], get_adaptive_inner_neighbor_);
		}
	}
	//=================================================================================================//
	SelfSurfaceContactRelation::
		SelfSurfaceContactRelation(RealBody &real_body)
		: BaseInnerRelation(real_body),
		  body_surface_layer_(real_body),
		  body_part_particles_(body_surface_layer_.body_part_particles_),
		  get_self_contact_neighbor_(real_body),
		  cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.getCellLinkedList())) {}
	//=================================================================================================//
	void SelfSurfaceContactRelation::resetNeighborhoodCurrentSize()
	{
		parallel_for(
			blocked_range<size_t>(0, body_part_particles_.size()),
			[&](const blocked_range<size_t> &r)
			{
				for (size_t num = r.begin(); num != r.end(); ++num)
				{
					size_t index_i = body_surface_layer_.getParticleIndex(num);
					inner_configuration_[index_i].current_size_ = 0;
				}
			},
			ap);
	}
	//=================================================================================================//
	void SelfSurfaceContactRelation::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = body_part_particles_.size();
		cell_linked_list_.searchNeighborsByParticles(
			body_surface_layer_, inner_configuration_,
			get_single_search_depth_, get_self_contact_neighbor_);
	}
	//=================================================================================================//
	void TreeInnerRelation::updateConfiguration()
	{
		generative_tree_.buildParticleConfiguration(inner_configuration_);
	}
	//=================================================================================================//
}