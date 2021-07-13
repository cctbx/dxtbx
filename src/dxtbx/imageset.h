/*
 * imageset.h
 *
 *  Copyright (C) 2017 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DXTBX_IMAGESET_H
#define DXTBX_IMAGESET_H

#include <map>

#include <boost/python.hpp>

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/sequence.h>
#include <dxtbx/format/image.h>
#include <dxtbx/error.h>
#include <dxtbx/masking/goniometer_shadow_masking.h>
#include <iostream>

namespace dxtbx {

using format::Image;
using format::ImageBuffer;
using format::ImageTile;
using masking::GoniometerShadowMasker;
using model::BeamBase;
using model::MonochromaticBeam;
using model::TOFBeam;
using model::Detector;
using model::Goniometer;
using model::Panel;
using model::Scan;
using model::TOFSequence;
using scitbx::rad_as_deg;
using scitbx::af::int2;

namespace detail {

  /**
   * Raise exception if we can't dereference pointer
   */
  template <typename T>
  T safe_dereference(boost::shared_ptr<T> ptr) {
    T *item = ptr.get();
    DXTBX_ASSERT(item != NULL);
    return *item;
  }

}  // namespace detail

/**
 * Hold an external lookup item
 */
template <typename T>
class ExternalLookupItem {
public:
  /** Construct the external lookup item */
  ExternalLookupItem() {}

  /**
   * Get the filename
   */
  std::string get_filename() const {
    return filename_;
  }

  /**
   * Set the filename
   * @param filename The filename
   */
  void set_filename(const std::string &filename) {
    filename_ = filename;
  }

  /**
   * Get the data
   */
  Image<T> get_data() const {
    return data_;
  }

  /**
   * Set the data
   * @param data - The data
   */
  void set_data(const Image<T> &data) {
    data_ = data;
  }

protected:
  std::string filename_;
  Image<T> data_;
};

/**
 * The external lookup group
 */
class ExternalLookup {
public:
  /**
   * @returns an internal reference to the mask
   */
  ExternalLookupItem<bool> &mask() {
    return mask_;
  }

  /**
   * @returns an internal reference to the gain
   */
  ExternalLookupItem<double> &gain() {
    return gain_;
  }

  /**
   * @returns an internal reference to the pedestal
   */
  ExternalLookupItem<double> &pedestal() {
    return pedestal_;
  }

  /**
   * @returns an internal reference to the dx map
   */
  ExternalLookupItem<double> &dx() {
    return dx_;
  }

  /**
   * @returns an internal reference to the dy map
   */
  ExternalLookupItem<double> &dy() {
    return dy_;
  }

protected:
  ExternalLookupItem<bool> mask_;
  ExternalLookupItem<double> gain_;
  ExternalLookupItem<double> pedestal_;
  ExternalLookupItem<double> dx_;
  ExternalLookupItem<double> dy_;
};

/**
 * A class to hold data about the imageset
 */
template<typename BeamT, typename SequenceT>
class ImageSetData {
public:
  typedef typename boost::shared_ptr<BeamT> beam_ptr;
  typedef boost::shared_ptr<Detector> detector_ptr;
  typedef boost::shared_ptr<Goniometer> goniometer_ptr;
  typedef boost::shared_ptr<SequenceT> sequence_ptr;
  typedef boost::shared_ptr<GoniometerShadowMasker> masker_ptr;

  ImageSetData() {}

  /**
   * Construct the imageset data object
   * @param reader The image reader
   * @param masker The image masker
   */
  ImageSetData(boost::python::object reader, masker_ptr masker)
      : reader_(reader),
        masker_(masker),
        beams_(boost::python::len(reader)),
        detectors_(boost::python::len(reader)),
        goniometers_(boost::python::len(reader)),
        sequences_(boost::python::len(reader)),
        reject_(boost::python::len(reader)) {}

  /**
   * @returns The reader object
   */
  boost::python::object reader() {
    return reader_;
  }

  /**
   * @returns The masker object
   */
  masker_ptr masker() {
    return masker_;
  }

  /**
   * @returns Does the imageset have a dynamic mask.
   */
  bool has_dynamic_mask() const {
    return masker_ != NULL;
  }

  /**
   * Read some image data
   * @param index The image index
   * @returns The image data
   */
  ImageBuffer get_data(std::size_t index) {
    // Create the return buffer
    ImageBuffer buffer;

    // Get the image data object
    boost::python::object data = reader_.attr("read")(index);

    // Get the class name
    std::string name =
      boost::python::extract<std::string>(data.attr("__class__").attr("__name__"))();

    // Extract the image buffer
    if (name == "tuple") {
      buffer = get_image_buffer_from_tuple(
        boost::python::extract<boost::python::tuple>(data)());
    } else {
      buffer = get_image_buffer_from_object(data);
    }
    return buffer;
  }

  /**
   * @returns Is the reader a single file reader
   */
  bool has_single_file_reader() const {
    return boost::python::extract<bool>(reader_.attr("is_single_file_reader")())();
  }

  /**
   * @returns The image path
   */
  std::string get_path(std::size_t index) const {
    return boost::python::extract<std::string>(reader_.attr("paths")()[index])();
  }

  /**
   * @returns The master image path
   */
  std::string get_master_path() const {
    return boost::python::extract<std::string>(reader_.attr("master_path")())();
  }

  /**
   * @returns The image identifier
   */
  std::string get_image_identifier(std::size_t index) const {
    return boost::python::extract<std::string>(reader_.attr("identifiers")()[index])();
  }

  /**
   * Mark the image for rejection
   * @param index The image index
   * @param reject True/False reject
   */
  void mark_for_rejection(std::size_t index, bool reject) {
    DXTBX_ASSERT(index < reject_.size());
    reject_[index] = reject;
  }

  /**
   * @param index The image index
   * @returns If the image is marked for rejection
   */
  bool is_marked_for_rejection(std::size_t index) const {
    DXTBX_ASSERT(index < reject_.size());
    return reject_[index];
  }

  /**
   * Get the array of booleans for image rejection
   */
  scitbx::af::shared<bool> get_reject_list() const {
    return reject_;
  }

  /**
   * Set the array of booleans for image rejection
   */
  void set_reject_list(const scitbx::af::const_ref<bool> &reject) {
    DXTBX_ASSERT(reject_.size() == reject.size());
    std::copy(reject.begin(), reject.end(), reject_.begin());
  }

  /**
   * Get a beam model
   * @param index The image index
   * @returns The beam model
   */
  beam_ptr get_beam(std::size_t index) const {
    DXTBX_ASSERT(index < beams_.size());
    return beams_[index];
  }

  /**
   * Set a beam model
   * @param index The image index
   * @param The beam model
   */
  void set_beam(const beam_ptr &beam, std::size_t index) {
    DXTBX_ASSERT(index < beams_.size());
    beams_[index] = beam;
  }

  /**
   * Get a detector model
   * @param index The image index
   * @returns The detector model
   */
  detector_ptr get_detector(std::size_t index) const {
    DXTBX_ASSERT(index < detectors_.size());
    return detectors_[index];
  }

  /**
   * Set a detector model
   * @param index The image index
   * @param The detector model
   */
  void set_detector(const detector_ptr &detector, std::size_t index) {
    DXTBX_ASSERT(index < detectors_.size());
    detectors_[index] = detector;
  }

  /**
   * Get a goniometer model
   * @param index The image index
   * @returns The goniometer model
   */
  goniometer_ptr get_goniometer(std::size_t index) const {
    DXTBX_ASSERT(index < goniometers_.size());
    return goniometers_[index];
  }

  /**
   * Set a goniometer model
   * @param index The image index
   * @param The goniometer model
   */
  void set_goniometer(const goniometer_ptr &goniometer, std::size_t index) {
    DXTBX_ASSERT(index < goniometers_.size());
    goniometers_[index] = goniometer;
  }

  /**
   * Get a sequence model
   * @param index The image index
   * @returns The sequence model
   */
  sequence_ptr get_sequence(std::size_t index) const {
    DXTBX_ASSERT(index < sequences_.size());
    return sequences_[index];
  }

  /**
   * Set a sequence model
   * @param index The image index
   * @param The sequence model
   */
  void set_sequence(const sequence_ptr &sequence, std::size_t index) {
    DXTBX_ASSERT(index < sequences_.size());
    sequences_[index] = sequence;
  }

  /**
   * @returns the number of images
   */
  std::size_t size() const {
    return boost::python::len(reader_);
  }

  /**
   * @returns The external lookup
   */
  ExternalLookup &external_lookup() {
    return external_lookup_;
  }

  /**
   * @returns the template
   */
  std::string get_template() const {
    return template_;
  }

  /**
   * @param x the template
   */
  void set_template(std::string x) {
    template_ = x;
  }

  /**
   * @returns the vendor
   */
  std::string get_vendor() const {
    return vendor_;
  }

  /**
   * @param x the vendor
   */
  void set_vendor(std::string x) {
    vendor_ = x;
  }

  /**
   * @returns the params
   */
  std::string get_params() const {
    return params_;
  }

  /**
   * @param x the params
   */
  void set_params(std::string x) {
    params_ = x;
  }

  /**
   * @returns the format
   */
  std::string get_format() const {
    return format_;
  }

  /**
   * @param x the format
   */
  void set_format(std::string x) {
    format_ = x;
  }

protected:
  ImageBuffer get_image_buffer_from_tuple(boost::python::tuple obj) {
    // Get the class name
    std::string name =
      boost::python::extract<std::string>(obj[0].attr("__class__").attr("__name__"))();

    // Read the image
    ImageBuffer buffer;
    if (name == "double") {
      buffer = ImageBuffer(get_image_from_tuple<double>(obj));
    } else if (name == "float") {
      buffer = ImageBuffer(get_image_from_tuple<float>(obj));
    } else if (name == "int") {
      buffer = ImageBuffer(get_image_from_tuple<int>(obj));
    } else {
      throw DXTBX_ERROR("Unknown type " + name);
    }
    return buffer;
  }

  ImageBuffer get_image_buffer_from_object(boost::python::object obj) {
    // Get the class name
    std::string name =
      boost::python::extract<std::string>(obj.attr("__class__").attr("__name__"))();

    // Read the image
    ImageBuffer buffer;
    if (name == "double") {
      buffer = ImageBuffer(Image<double>(get_image_tile_from_object<double>(obj)));
    } else if (name == "float") {
      buffer = ImageBuffer(Image<float>(get_image_tile_from_object<float>(obj)));
    } else if (name == "int") {
      buffer = ImageBuffer(Image<int>(get_image_tile_from_object<int>(obj)));
    } else {
      throw DXTBX_ERROR("Unknown type " + name);
    }
    return buffer;
  }

  template <typename T>
  Image<T> get_image_from_tuple(boost::python::tuple obj) {
    Image<T> image;
    for (std::size_t i = 0; i < boost::python::len(obj); ++i) {
      image.push_back(get_image_tile_from_object<T>(obj[i]));
    }
    return image;
  }

  template <typename T>
  ImageTile<T> get_image_tile_from_object(boost::python::object obj) {
    typedef typename scitbx::af::flex<T>::type flex_type;

    // Extract the data
    flex_type a = boost::python::extract<flex_type>(obj)();

    // Return the image tile
    return ImageTile<T>(scitbx::af::versa<T, scitbx::af::c_grid<2> >(
      a.handle(), scitbx::af::c_grid<2>(a.accessor())));
  }

  boost::python::object reader_;
  boost::shared_ptr<GoniometerShadowMasker> masker_;
  scitbx::af::shared<beam_ptr> beams_;
  scitbx::af::shared<detector_ptr> detectors_;
  scitbx::af::shared<goniometer_ptr> goniometers_;
  scitbx::af::shared<sequence_ptr> sequences_;
  scitbx::af::shared<bool> reject_;
  ExternalLookup external_lookup_;

  std::string template_;
  std::string vendor_;
  std::string params_;
  std::string format_;
};

/**
 * A class to represent an imageset
 */
template<class BeamT, class SequenceT>
class ImageSet {
public:
  typedef typename ImageSetData<BeamT, SequenceT>::beam_ptr beam_ptr;
  typedef typename ImageSetData<BeamT, SequenceT>::detector_ptr detector_ptr;
  typedef typename ImageSetData<BeamT, SequenceT>::goniometer_ptr goniometer_ptr;
  typedef typename ImageSetData<BeamT, SequenceT>::sequence_ptr sequence_ptr;

  /**
   * Cache an image
   */
  template <class T>
  class DataCache {
  public:
    T image;
    int index;

    DataCache() : index(-1) {}

    bool image_in_cache(std::size_t idx){
      return index == idx; 
    }

  };

  /**
   * Default constructor throws an exception.
   * This only here so overloaded functions that
   * never return can be compiled without warnings.
   */
  ImageSet() {
    throw DXTBX_ERROR("ImageSet needs imageset data");
  }

  /**
   * Destructor
   */
  virtual ~ImageSet() {}

  /**
   * Construct the imageset
   * @param data The imageset data
   */
  ImageSet(const ImageSetData<BeamT, SequenceT> &data) : data_(data), indices_(data.size()) {
    // Check number of images
    if (data.size() == 0) {
      throw DXTBX_ERROR("No images specified in ImageSetData");
    }

    // Set indices
    for (std::size_t i = 0; i < indices_.size(); ++i) {
      indices_[i] = i;
    }
  }

  /**
   * Construct the imageset
   * @param data The imageset data
   * @param indices The image indices
   */
  ImageSet(const ImageSetData<BeamT, SequenceT> &data, const scitbx::af::const_ref<std::size_t> &indices)
      : data_(data), indices_(indices.begin(), indices.end()) {
    // Check number of images
    if (data.size() == 0) {
      throw DXTBX_ERROR("No images specified in ImageSetData");
    }

    // Set the indices
    if (indices.size() == 0) {
      indices_.resize(data.size());
      for (std::size_t i = 0; i < indices_.size(); ++i) {
        indices_[i] = i;
      }
    } else {
      // Check indices
      if (scitbx::af::max(indices) >= data.size()) {
        throw DXTBX_ERROR("Indices are not consistent with # of images");
      }
    }
  }

  /**
   * @returns the imageset data
   */
  ImageSetData<BeamT, SequenceT> data() const {
    return data_;
  }

  /**
   * @returns The image indices
   */
  scitbx::af::shared<std::size_t> indices() const {
    return indices_;
  }

  /**
   * @returns The number of images
   */
  std::size_t size() const {
    return indices_.size();
  }

  /**
   * @returns The external lookup
   */
  ExternalLookup &external_lookup() {
    return data_.external_lookup();
  }

  /**
   * Get the raw image data
   * @param index The image index
   * @returns The raw image data
   */
  ImageBuffer get_raw_data(std::size_t index) {
    DXTBX_ASSERT(index < indices_.size());
    if (data_cache_.index == index) {
      return data_cache_.image;
    }
    ImageBuffer image = data_.get_data(indices_[index]);
    data_cache_.index = index;
    data_cache_.image = image;
    return image;
  }

  /**
   * Get the corrected data array (raw - pedestal) / gain
   * @param index The image index
   * @returns The corrected data array
   */
  virtual Image<double> get_corrected_data(std::size_t index) {
    typedef scitbx::af::versa<double, scitbx::af::c_grid<2> > array_type;
    typedef scitbx::af::const_ref<double, scitbx::af::c_grid<2> > const_ref_type;

    // Get the multi-tile data, gain and pedestal
    DXTBX_ASSERT(index < indices_.size());
    Image<double> data = get_raw_data_as_double(index);
    Image<double> gain = get_gain(index);
    Image<double> dark = get_pedestal(index);
    DXTBX_ASSERT(gain.n_tiles() == 0 || data.n_tiles() == gain.n_tiles());
    DXTBX_ASSERT(dark.n_tiles() == 0 || data.n_tiles() == dark.n_tiles());

    // Loop through tiles
    Image<double> result;
    for (std::size_t i = 0; i < data.n_tiles(); ++i) {
      // Get the data
      const_ref_type r = data.tile(i).data().const_ref();

      // Get the gain and dark
      const_ref_type g = gain.n_tiles() > 0
                           ? gain.tile(i).data().const_ref()
                           : const_ref_type(NULL, scitbx::af::c_grid<2>(0, 0));
      const_ref_type p = dark.n_tiles() > 0
                           ? dark.tile(i).data().const_ref()
                           : const_ref_type(NULL, scitbx::af::c_grid<2>(0, 0));

      // Check gain and dark sizes
      DXTBX_ASSERT(g.size() == 0 || r.accessor().all_eq(g.accessor()));
      DXTBX_ASSERT(p.size() == 0 || r.accessor().all_eq(p.accessor()));

      if (p.size() == 0 && g.size() == 0) {
        // Nothing to apply, save the copy
        result.push_back(ImageTile<double>(data.tile(i).data()));
      } else {
        // Create the result array
        array_type c(r.accessor(),
                     scitbx::af::init_functor_null<array_type::value_type>());

        // Copy the data values
        std::uninitialized_copy(r.begin(), r.end(), c.begin());

        // Apply dark
        if (p.size() > 0) {
          for (std::size_t j = 0; j < r.size(); ++j) {
            c[j] = c[j] - p[j];
          }
        }

        // Apply gain
        if (g.size() > 0) {
          for (std::size_t j = 0; j < r.size(); ++j) {
            DXTBX_ASSERT(g[j] > 0);
            c[j] = c[j] / g[j];
          }
        }

        // Add the image tile
        result.push_back(ImageTile<double>(c));
      }
    }

    // Return the result
    return result;
  }

  /**
   * Get the detector gain map. Either take this from the external gain map or
   * try to construct from the detector gain.
   * @param index The image index
   * @returns The gain
   */
  Image<double> get_gain(std::size_t index) {
    // If the external lookup is empty
    DXTBX_ASSERT(index < indices_.size());
    if (external_lookup().gain().get_data().empty()) {
      // Get the detector
      Detector detector = detail::safe_dereference(get_detector_for_image(index));

      // Compute the gain for each panel
      bool use_detector_gain = true;
      bool need_gain_map = false;
      std::vector<double> gain(detector.size(), 0);
      for (std::size_t i = 0; i < detector.size(); ++i) {
        gain[i] = detector[i].get_gain();
        if (gain[i] <= 0) {
          use_detector_gain = false;
          break;
        } else if (std::abs(gain[i] - 1.0) > 1e-7) {
          need_gain_map = true;
        }
      }

      // If using the gain from the panel, construct a gain map
      if (use_detector_gain && need_gain_map) {
        Image<double> result;
        for (std::size_t i = 0; i < detector.size(); ++i) {
          std::size_t xsize = detector[i].get_image_size()[0];
          std::size_t ysize = detector[i].get_image_size()[1];
          scitbx::af::c_grid<2> grid(ysize, xsize);
          scitbx::af::versa<double, scitbx::af::c_grid<2> > data(grid, gain[i]);
          result.push_back(ImageTile<double>(data));
        }
        return result;
      }
    }
    return external_lookup().gain().get_data();
  }

  /**
   * Get the pedestal
   * @param index The image index
   * @returns The pedestal image
   */
  Image<double> get_pedestal(std::size_t index) {
    //
    // If the external lookup is empty
    DXTBX_ASSERT(index < indices_.size());
    if (external_lookup().pedestal().get_data().empty()) {
      // Get the detector
      Detector detector = detail::safe_dereference(get_detector_for_image(index));

      // Compute the pedestal for each panel
      bool use_detector_pedestal = false;
      std::vector<double> pedestal(detector.size(), 0);
      for (std::size_t i = 0; i < detector.size(); ++i) {
        pedestal[i] = detector[i].get_pedestal();
        if (std::abs(pedestal[i]) > 1e-7) {
          use_detector_pedestal = true;
        }
      }

      // If using the pedestal from the panel, construct a pedestal map
      if (use_detector_pedestal) {
        Image<double> result;
        for (std::size_t i = 0; i < detector.size(); ++i) {
          std::size_t xsize = detector[i].get_image_size()[0];
          std::size_t ysize = detector[i].get_image_size()[1];
          scitbx::af::c_grid<2> grid(ysize, xsize);
          scitbx::af::versa<double, scitbx::af::c_grid<2> > data(grid, pedestal[i]);
          result.push_back(ImageTile<double>(data));
        }
        return result;
      }
    }
    return external_lookup().pedestal().get_data();
  }

  /**
   * @returns Does the imageset have a dynamic mask.
   */
  bool has_dynamic_mask() const {
    return data_.has_dynamic_mask();
  }

  /**
   * Get an empty mask
   * @param index The image index
   * @returns The mask
   */
  Image<bool> get_empty_mask() const {
    Detector detector = detail::safe_dereference(get_detector_for_image(0));
    Image<bool> mask;
    for (std::size_t i = 0; i < detector.size(); ++i) {
      std::size_t xsize = detector[i].get_image_size()[0];
      std::size_t ysize = detector[i].get_image_size()[1];
      mask.push_back(ImageTile<bool>(scitbx::af::versa<bool, scitbx::af::c_grid<2> >(
        scitbx::af::c_grid<2>(ysize, xsize), true)));
    }
    return mask;
  }

  /**
   * Get the untrusted rectangle mask
   * @param mask The mask to write into
   * @returns The mask
   */
  Image<bool> get_untrusted_rectangle_mask(Image<bool> mask) const {
    Detector detector = detail::safe_dereference(get_detector_for_image(0));
    DXTBX_ASSERT(mask.n_tiles() == detector.size());
    for (std::size_t i = 0; i < detector.size(); ++i) {
      detector[i].apply_untrusted_rectangle_mask(mask.tile(i).data().ref());
    }
    return mask;
  }

  /**
   * Apply the external mask
   * @param mask The input mask
   * @returns The external mask
   */
  Image<bool> get_external_mask(Image<bool> mask) {
    Image<bool> external_mask = external_lookup().mask().get_data();
    if (!external_mask.empty()) {
      DXTBX_ASSERT(external_mask.n_tiles() == mask.n_tiles());
      for (std::size_t i = 0; i < mask.n_tiles(); ++i) {
        scitbx::af::ref<bool, scitbx::af::c_grid<2> > m1 = mask.tile(i).data().ref();
        scitbx::af::const_ref<bool, scitbx::af::c_grid<2> > m2 =
          external_mask.tile(i).data().const_ref();
        DXTBX_ASSERT(m1.accessor().all_eq(m2.accessor()));
        for (std::size_t j = 0; j < m1.size(); ++j) {
          m1[j] = m1[j] && m2[j];
        }
      }
    }
    return mask;
  }

  /**
   * Get the static mask common to all images
   * @param mask The input mask
   * @returns The mask
   */
  Image<bool> get_static_mask(Image<bool> mask) {
    return get_untrusted_rectangle_mask(
      get_external_mask(mask.empty() ? get_empty_mask() : mask));
  }

  /**
   * Get the static mask common to all images
   * @returns The mask
   */
  Image<bool> get_static_mask() {
    return get_static_mask(Image<bool>());
  }

  /**
   * Get the trusted range mask for the index
   * @param mask The mask to write into
   * @param index The image index
   * @returns The mask
   */
  Image<bool> get_trusted_range_mask(Image<bool> mask, std::size_t index) {
    Detector detector = detail::safe_dereference(get_detector_for_image(index));
    Image<double> data = get_raw_data_as_double(index);
    DXTBX_ASSERT(mask.n_tiles() == data.n_tiles());
    DXTBX_ASSERT(data.n_tiles() == detector.size());
    for (std::size_t i = 0; i < detector.size(); ++i) {
      detector[i].apply_trusted_range_mask(data.tile(i).data().const_ref(),
                                           mask.tile(i).data().ref());
    }
    return mask;
  }

  /**
   * Get the dynamic mask for the requested image
   * @param index The image index
   * @returns The image mask
   */
  virtual Image<bool> get_dynamic_mask(std::size_t index) {
    return get_trusted_range_mask(get_static_mask(), index);
  }

  /**
   * Compute the mask
   * @param index The image index
   * @returns The image mask
   */
  Image<bool> get_mask(std::size_t index) {
    return get_dynamic_mask(index);
  }

  /**
   * @param index The image index
   * @returns the beam at index
   */
  virtual beam_ptr get_beam_for_image(std::size_t index = 0) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.get_beam(indices_[index]);
  }

  /**
   * @param index The image index
   * @returns the detector at index
   */
  virtual detector_ptr get_detector_for_image(std::size_t index = 0) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.get_detector(indices_[index]);
  }

  /**
   * @param index The image index
   * @returns the goniometer at index
   */
  virtual goniometer_ptr get_goniometer_for_image(std::size_t index = 0) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.get_goniometer(indices_[index]);
  }

  /**
   * @param index The image index
   * @returns the sequence at index
   */
  virtual sequence_ptr get_sequence_for_image(std::size_t index = 0) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.get_sequence(indices_[index]);
  }

  /**
   * Set the beam model
   * @param index The image index
   * @param beam The beam model
   */
  virtual void set_beam_for_image(const beam_ptr &beam, std::size_t index = 0) {
    DXTBX_ASSERT(index < indices_.size());
    data_.set_beam(beam, indices_[index]);
  }

  /**
   * Set the detector model
   * @param index The image index
   * @param detector The detector model
   */
  virtual void set_detector_for_image(const detector_ptr &detector,
                                      std::size_t index = 0) {
    DXTBX_ASSERT(index < indices_.size());
    data_.set_detector(detector, indices_[index]);
  }

  /**
   * Set the goniometer model
   * @param index The image index
   * @param goniometer The goniometer model
   */
  virtual void set_goniometer_for_image(const goniometer_ptr &goniometer,
                                        std::size_t index = 0) {
    DXTBX_ASSERT(index < indices_.size());
    data_.set_goniometer(goniometer, indices_[index]);
  }

  /**
   * Set the sequence model
   * @param index The image index
   * @param sequence The sequence model
   */
  virtual void set_sequence_for_image(const sequence_ptr &sequence, std::size_t index = 0) {
    DXTBX_ASSERT(sequence == NULL || sequence->get_num_images() == 1);
    DXTBX_ASSERT(index < indices_.size());
    data_.set_sequence(sequence, indices_[index]);
  }

  /**
   * Get the image path
   * @param index The image index
   * @returns The image path
   */
  std::string get_path(std::size_t index) const {
    DXTBX_ASSERT(index < indices_.size());
    if (data_.has_single_file_reader()) {
      return data_.get_master_path();
    }
    return data_.get_path(indices_[index]);
  }

  /**
   * Get the image identifier
   * @param index The image index
   * @returns The image identifier
   */
  std::string get_image_identifier(std::size_t index) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.get_image_identifier(indices_[index]);
  }

  /**
   * Mark the image for rejection
   * @param index The image index
   * @param reject True/False reject
   */
  void mark_for_rejection(std::size_t index, bool reject) {
    DXTBX_ASSERT(index < indices_.size());
    data_.mark_for_rejection(indices_[index], reject);
  }

  /**
   * @param index The image index
   * @returns If the image is marked for rejection
   */
  bool is_marked_for_rejection(std::size_t index) const {
    DXTBX_ASSERT(index < indices_.size());
    return data_.is_marked_for_rejection(indices_[index]);
  }

  /**
   * @returns The imageset itself
   */
  virtual ImageSet as_imageset() const {
    return *this;
  }

  /**
   * @returns The complete set
   */
  virtual ImageSet complete_set() const {
    return ImageSet(data_);
  }

  /**
   * @param first The first slice index
   * @param last The last slice index
   * @returns The partial set
   */
  ImageSet partial_set(std::size_t first, std::size_t last) const {
    DXTBX_ASSERT(last > first);
    return ImageSet(data_,
                    scitbx::af::const_ref<std::size_t>(&indices_[first], last - first));
  }

  /**
   * Compare the imageset with another
   * @param other The other imageset
   * @returns Are the imagesets the same
   */
  bool operator==(const ImageSet &other) const {
    if (size() == other.size()) {
      for (std::size_t i = 0; i < size(); ++i) {
        if (get_path(i) != other.get_path(i)) {
          return false;
        }
        if (indices()[i] != other.indices()[i]) {
          return false;
        }
      }
      return true;
    }
    return false;
  }

  /**
   * Compare the imageset with another
   */
  bool operator!=(const ImageSet &other) const {
    return !(*this == other);
  }

  /**
   * Clear the imageset cache. Useful for when many imagesets are
   * held in memory. After reading an image the cache can be
   * manually cleared before moving onto the next imageset.
   */
  void clear_cache() {
    data_cache_ = DataCache<ImageBuffer>();
    double_raw_data_cache_ = DataCache<Image<double> >();
  }

protected:
  ImageSetData<BeamT, SequenceT> data_;
  scitbx::af::shared<std::size_t> indices_;
  DataCache<ImageBuffer> data_cache_;
  DataCache<Image<double> > double_raw_data_cache_;

  Image<double> get_raw_data_as_double(std::size_t index) {
    DXTBX_ASSERT(index < indices_.size());
    if (double_raw_data_cache_.image_in_cache(index)) {
      return double_raw_data_cache_.image;
    }
    Image<double> image = get_raw_data(index).as_double();
    double_raw_data_cache_.index = index;
    double_raw_data_cache_.image = image;
    return image;
  }

};

/**
 * Class to represent a grid of images
 */
class ImageGrid : public ImageSet<MonochromaticBeam, Scan> {
public:
  /**
   * Construct the grid
   * @param data The imageset data
   * @param grid_size The size of the grid
   */
  ImageGrid(const ImageSetData<MonochromaticBeam, Scan> &data, int2 grid_size)
      : ImageSet<MonochromaticBeam, Scan>(data), grid_size_(grid_size) {
    DXTBX_ASSERT(grid_size.all_gt(0));
    DXTBX_ASSERT(grid_size[0] * grid_size[1] == size());
  }

  /**
   * Construct the grid
   * @param data The imageset data
   * @param indices The imageset indices
   * @param grid_size The size of the grid
   */
  ImageGrid(const ImageSetData<MonochromaticBeam, Scan> &data,
            const scitbx::af::const_ref<std::size_t> &indices,
            int2 grid_size)
      : ImageSet<MonochromaticBeam, Scan>(data, indices), grid_size_(grid_size) {
    DXTBX_ASSERT(grid_size.all_gt(0));
    DXTBX_ASSERT(grid_size[0] * grid_size[1] == size());
  }

  /**
   * Destructor
   */
  virtual ~ImageGrid() {}

  /**
   * @returns The grid size
   */
  int2 get_grid_size() const {
    return grid_size_;
  }

  /**
   * Construct the image grid from the imageset
   * @param imageset The imageset
   * @param grid_size The grid_size
   * @returns the image grid
   */
  static ImageGrid from_imageset(const ImageSet<MonochromaticBeam, Scan> &imageset, int2 grid_size) {
    ImageGrid result(imageset.data(), imageset.indices().const_ref(), grid_size);
    return result;
  }

  /**
   * Convert the grid to an imageset
   * @returns An imageset
   */
  ImageSet<MonochromaticBeam, Scan> as_imageset() const {
    ImageSet<MonochromaticBeam, Scan> result(data_, indices_.const_ref());
    return result;
  }

  /**
   * Get the complete set
   * @returns The complete sequence
   */
  ImageSet<MonochromaticBeam, Scan> complete_set() const {
    throw DXTBX_ERROR("Cannot get complete set from image grid");
    return ImageSet();
  }

  /**
   * Get a partial set
   * @param first The first index
   * @param last The last index
   * @returns The partial sequence
   */
  ImageSet<MonochromaticBeam, Scan> partial_set(std::size_t first, std::size_t last) const {
    throw DXTBX_ERROR("Cannot get partial set from image grid");
    return ImageSet<MonochromaticBeam, Scan>();
  }

protected:
  int2 grid_size_;
};


/**
 * A class to represent a sequence of data
 */
template<class BeamT, class SequenceT>
class ImageSequence : public ImageSet<BeamT, SequenceT> {
public:
  typedef typename ImageSet<BeamT, SequenceT>::beam_ptr beam_ptr;
  typedef typename ImageSet<BeamT, SequenceT>::detector_ptr detector_ptr;
  typedef typename ImageSet<BeamT, SequenceT>::goniometer_ptr goniometer_ptr;
  typedef typename ImageSet<BeamT, SequenceT>::sequence_ptr sequence_ptr;
  /**
   * Construct the sequence
   * @param data The imageset data
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  ImageSequence(const ImageSetData<BeamT, SequenceT> &data,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSet<BeamT , SequenceT>(data),
        beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        sequence_(sequence) {
    // Check the sequence is the same length and number of images
    DXTBX_ASSERT(sequence.get() != NULL);
    if (data.size() > 1) {
      if (data.size() != sequence->get_num_images()) {
        throw DXTBX_ERROR("Sequence size is not compatible with number of images");
      }
    }

    // Set the models for each image
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_beam_for_image(beam_, i);
      ImageSet<BeamT, SequenceT>::set_detector_for_image(detector_, i);
      ImageSet<BeamT, SequenceT>::set_goniometer_for_image(goniometer_, i);
      ImageSet<BeamT, SequenceT>::set_sequence_for_image(sequence_ptr(new SequenceT((*sequence)[i])), i);
    }
  }

  /**
   * Construct the sequence
   * @param data The imageset data
   * @param indices The image indices
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  ImageSequence(const ImageSetData<BeamT, SequenceT> &data,
                const scitbx::af::const_ref<std::size_t> &indices,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSet<BeamT, SequenceT>(data, indices),
        beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        sequence_(sequence) {
    // Check the sequence is the same length as number of indices
    DXTBX_ASSERT(sequence.get() != NULL);

    // Check indices are sequential
    if (indices.size() > 0) {
      if (indices.size() != sequence->get_num_images()) {
        throw DXTBX_ERROR("Sequence size is not compatible with number of images");
      }
      for (std::size_t i = 1; i < indices.size(); ++i) {
        DXTBX_ASSERT(indices[i] == indices[i - 1] + 1);
      }
    }

    // Set the models for each image
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_beam_for_image(beam_, i);
      ImageSet<BeamT, SequenceT>::set_detector_for_image(detector_, i);
      ImageSet<BeamT, SequenceT>::set_goniometer_for_image(goniometer_, i);
      ImageSet<BeamT, SequenceT>::set_sequence_for_image(sequence_ptr(new SequenceT((*sequence)[i])), i);
    }
  }

  /**
   * Destructor
   */
  virtual ~ImageSequence() {}


  /**
   * @returns the array range
   */
  int2 get_array_range() const {
    DXTBX_ASSERT(sequence_ != NULL);
    return sequence_->get_array_range();
  }

  /**
   * @returns the beam model
   */
  beam_ptr get_beam() const {
    return beam_;
  }

  /**
   * @returns the detector model
   */
  detector_ptr get_detector() const {
    return detector_;
  }

  /**
   * @returns the goniometer model
   */
  goniometer_ptr get_goniometer() const {
    return goniometer_;
  }

  /**
   * @returns the sequence model
   */
  sequence_ptr get_sequence() const {
    return sequence_;
  }

  /**
   * Set the beam model
   * @param beam The beam model
   */
  void set_beam(const beam_ptr &beam) {
    beam_ = beam;
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_beam_for_image(beam_, i);
    }
  }

  /**
   * Set the detector model
   * @param detector The detector model
   */
  void set_detector(const detector_ptr &detector) {
    detector_ = detector;
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_detector_for_image(detector_, i);
    }
  }

  /**
   * Set the goniometer model
   * @param goniometer The goniometer model
   */
  void set_goniometer(const goniometer_ptr &goniometer) {
    goniometer_ = goniometer;
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_goniometer_for_image(goniometer_, i);
    }
  }

  /**
   * Set the sequence model
   * @param sequence The sequence model
   */
  void set_sequence(const sequence_ptr &sequence) {
    DXTBX_ASSERT(sequence.get() != NULL);
    if (sequence->get_num_images() != ImageSet<BeamT, SequenceT>::size()) {
      DXTBX_ASSERT(sequence_ != NULL);
      int i0 = sequence->get_array_range()[0];
      int i1 = sequence->get_array_range()[1];
      int j0 = sequence_->get_array_range()[0];
      DXTBX_ASSERT(i0 >= j0);
      DXTBX_ASSERT(i1 > i0);
      std::size_t n = i1 - i0;
      int k0 = i0 - j0;
      DXTBX_ASSERT(k0 >= 0);
      std::size_t index0 = ImageSet<BeamT, SequenceT>::indices_[0];
      ImageSet<BeamT, SequenceT>::indices_.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        ImageSet<BeamT, SequenceT>::indices_[i] = index0 + i;
      }
    }
    DXTBX_ASSERT((sequence->get_num_images() == ImageSet<BeamT, SequenceT>::size()));
    sequence_ = sequence;
    for (std::size_t i = 0; i < ImageSet<BeamT, SequenceT>::size(); ++i) {
      ImageSet<BeamT, SequenceT>::set_sequence_for_image(sequence_ptr(new SequenceT((*sequence)[i])), i);
    }
  }

  /**
   * Override per-image model
   */
  void set_beam_for_image(const beam_ptr &beam, std::size_t index) {
    throw DXTBX_ERROR("Cannot set per-image model in sequence");
  }

  /**
   * Override per-image model
   */
  void set_detector_for_image(const detector_ptr &detector, std::size_t index) {
    throw DXTBX_ERROR("Cannot set per-image model in sequence");
  }

  /**
   * Override per-image model
   */
  void set_goniometer_for_image(const goniometer_ptr &goniometer, std::size_t index) {
    throw DXTBX_ERROR("Cannot set per-image model in sequence");
  }

  /**
   * Override per-image model
   */
  void set_sequence_for_image(const sequence_ptr &sequence, std::size_t index) {
    throw DXTBX_ERROR("Cannot set per-image model in sequence");
  }

  /**
   * Convert the sequence to an imageset
   * @returns An imageset
   */
  ImageSet<BeamT, SequenceT> as_imageset() const {
    ImageSet<BeamT, SequenceT> result(ImageSet<BeamT, SequenceT>::data_, 
                                    ImageSet<BeamT, SequenceT>::indices_.const_ref());
    return result;
  }

  /**
   * Get the complete seta
   * @returns The complete sequence
   */
  ImageSet<BeamT, SequenceT> complete_set() const {
    throw DXTBX_ERROR("Cannot get complete set from image sequence");
    return ImageSet<BeamT, SequenceT>();
  }

  /**
   * Get a partial set
   * @param first The first index
   * @param last The last index
   * @returns The partial sequence
   */
  ImageSet<BeamT, SequenceT> partial_set(std::size_t first, std::size_t last) const {
    throw DXTBX_ERROR("Cannot get partial set from image sequence");
    return ImageSet<BeamT, SequenceT>();
  }

  /**
   * Get the complete seta
   * @returns The complete sequence
   */
  ImageSequence<BeamT, SequenceT> complete_sequence() const {
    // Compute sequence
    SequenceT sequence = detail::safe_dereference(ImageSet<BeamT, SequenceT>::data_.get_sequence(0));
    for (std::size_t i = 1; i < ImageSet<BeamT, SequenceT>::data_.size(); ++i) {
      sequence += detail::safe_dereference(ImageSet<BeamT, SequenceT>::data_.get_sequence(i));
    }
   
    // Construct a sequence
    ImageSequence<BeamT, SequenceT> result(
      ImageSet<BeamT, SequenceT>::data_, get_beam(), get_detector(), get_goniometer(), sequence_ptr(new SequenceT(sequence)));

    // Return the sequence
    return result;
  }

  /**
   * Get a partial set
   * @param first The first index
   * @param last The last index
   * @returns The partial sequence
   */
  ImageSequence<BeamT, SequenceT> partial_sequence(std::size_t first, std::size_t last) const {
    // Check slice indices
    DXTBX_ASSERT(last > first);

    // Construct a partial sequence
    SequenceT sequence = detail::safe_dereference(ImageSet<BeamT, SequenceT>::get_sequence_for_image(first));
    for (std::size_t i = first + 1; i < last; ++i) {
      sequence += detail::safe_dereference(ImageSet<BeamT, SequenceT>::get_sequence_for_image(i));
    }

    // Construct the partial indices
    scitbx::af::const_ref<std::size_t> indices(&ImageSet<BeamT, SequenceT>::indices_[first], last - first);

    // Construct the partial sequence
    ImageSequence<BeamT, SequenceT> result(ImageSet<BeamT, SequenceT>::data_,
                         indices,
                         get_beam(),
                         get_detector(),
                         get_goniometer(),
                         sequence_ptr(new SequenceT(sequence)));

    // Return the sequence
    return result;
  }

protected:
  beam_ptr beam_;
  detector_ptr detector_;
  goniometer_ptr goniometer_;
  sequence_ptr sequence_;
};

/**
 * A class to represent images with a ToF dimension
 */
class TOFImageSequence : public ImageSequence<TOFBeam, TOFSequence>{

public:

  /**
   * @param data The imageset data
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  TOFImageSequence(const ImageSetData<TOFBeam, TOFSequence> &data,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSequence<TOFBeam, TOFSequence>(data, beam, detector, goniometer, sequence){}

  /**
   * @param data The imageset data
   * @param indices The image indices
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  TOFImageSequence(const ImageSetData<TOFBeam, TOFSequence> &data,
                const scitbx::af::const_ref<std::size_t> &indices,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSequence<TOFBeam, TOFSequence>(data, indices, beam, detector, goniometer, sequence){}

  virtual ~TOFImageSequence() {}
  
  scitbx::af::shared<double> tof_in_seconds() const{
    return sequence_->get_tof_in_seconds();
  }

};

class RotImageSequence : public ImageSequence<MonochromaticBeam, Scan>{
public:

  /**
   * @param data The imageset data
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  RotImageSequence(const ImageSetData<MonochromaticBeam, Scan> &data,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSequence<MonochromaticBeam, Scan>(data, beam, detector, goniometer, sequence){}

  /**
   * @param data The imageset data
   * @param indices The image indices
   * @param beam The beam model
   * @param detector The detector model
   * @param goniometer The gonioeter model
   * @param sequence The sequence model
   */
  RotImageSequence(const ImageSetData<MonochromaticBeam, Scan> &data,
                const scitbx::af::const_ref<std::size_t> &indices,
                const beam_ptr &beam,
                const detector_ptr &detector,
                const goniometer_ptr &goniometer,
                const sequence_ptr &sequence)
      : ImageSequence<MonochromaticBeam, Scan>(data, indices, beam, detector, goniometer, sequence){}

  virtual ~RotImageSequence() {}
  
  /**
   * Get the dynamic mask for the requested image
   * @param index The image index
   * @returns The image mask
   */
  virtual Image<bool> get_dynamic_mask(std::size_t index) {
    // Get the masker
    ImageSetData<MonochromaticBeam, Scan>::masker_ptr masker = ImageSet<MonochromaticBeam, Scan>::data_.masker();

    // Create return buffer
    Image<bool> dyn_mask;

    // Get the image data object
    if (masker != NULL) {
      DXTBX_ASSERT(sequence_ != NULL);
      DXTBX_ASSERT(detector_ != NULL);
      double scan_angle = rad_as_deg(
        sequence_->get_angle_from_image_index(index + sequence_->get_image_range()[0]));
      dyn_mask = masker->get_mask(*detector_, scan_angle);
    }

    // Return the dynamic mask
    return get_trusted_range_mask(get_static_mask(dyn_mask), index);
  }
  

};

}  // namespace dxtbx

#endif  // DXTBX_IMAGESET_H
