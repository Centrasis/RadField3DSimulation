#pragma once
#include <memory>
#include <string>
#include <radfiled3d/radiation_field.hpp>
#include <radfiled3d/storage/types.hpp>


namespace RadiationSimulation {
	/** Appends a field simulated by RadField3D to a RadField3D field file, so that the file holds the result of all joined
	* runs as if they had been one run.
	*
	* Both fields are normalized per primary particle, so every channel of `field` is combined with the file's channel of
	* the same name, weighting each run by its share r of the primaries:
	* - flux, angular_flux: (1 - r) * existing + r * new
	* - spectrum: per voxel, the spectra mixed by each run's flux contribution and renormalized to sum 1
	* - error: sqrt((w_existing e_existing)^2 + (w_new e_new)^2) / (w_existing + w_new) with the same per-voxel flux weights
	* The primary particle counts and simulation durations add up. Two runs with the same random seed are refused (their
	* random numbers would repeat each other); the combined file records the random seed 0, as it stems from several runs. Channels of the file that `field` does not contain
	* (e.g. the voxelized geometry) are kept unchanged. All arithmetic runs in double precision; a result that does not
	* fit the stored type, an unknown layer, or a mismatching grid, layer layout or simulation setup raise an exception and
	* leave the file untouched. If the file does not exist, `field` is stored as it is.
	*
	* The file is locked while it is updated and replaced atomically through a temporary file next to it.
	* @param field The field to append, normalized per primary particle.
	* @param metadata The metadata of `field`.
	* @param file The field file to append to.
	*/
	void append_radiation_field(std::shared_ptr<radfiled3d::CartesianRadiationField> field, std::shared_ptr<radfiled3d::storage::v1::RadiationFieldMetadata> metadata, const std::string& file);

	/** Stores a field so that `file` is never seen incomplete: it is written to a temporary file next to `file`, which then
	* replaces `file` in a single rename.
	* @param field The field to store.
	* @param metadata The metadata of `field`.
	* @param file The file to create or overwrite.
	*/
	void store_radiation_field_atomically(std::shared_ptr<radfiled3d::IRadiationField> field, std::shared_ptr<radfiled3d::storage::RadiationFieldMetadata> metadata, const std::string& file);

	/** A token unique to this process and call, for names of files only this process writes. It is built from the
	* process id, the clock and a counter, not from a random engine, so naming files never consumes random numbers of
	* the simulation.
	*/
	std::string unique_file_token();
}
