#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/shared_ptr.hpp>
#include "Python.h"

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <vector>
#include <dxtbx/imageset.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/sequence.h>
#include <dxtbx/model/pixel_to_millimeter.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace boost_python {

  using model::OffsetParallaxCorrectedPxMmStrategy;
  using model::OffsetPxMmStrategy;

  namespace detail {

    /// Returns a byte-data-appropriate object for an std::string
    boost::python::object bytes_from_std_string(const std::string &data) {
      // Boost::python appears to have no native way to do this conversion
      return boost::python::object(
        boost::python::handle<>(PyBytes_FromStringAndSize(data.data(), data.size())));
    }

    /// Pickle a python object to a string
    std::string pickle_dumps(boost::python::object x) {
      return boost::python::extract<std::string>(
        boost::python::import("pickle").attr("dumps")(x));
    }

    /// Unpickle a python object from a string
    boost::python::object pickle_loads(std::string x) {
      if (x == "") {
        return boost::python::object();
      }

      return boost::python::import("pickle").attr("loads")(bytes_from_std_string(x));
    }
  }  // namespace detail

  template<typename BeamT, typename SequenceT>
  typename ImageSetData<BeamT, SequenceT>::masker_ptr make_masker_pointer(boost::python::object masker) {
    if (masker == boost::python::object()) {
      return typename ImageSetData<BeamT, SequenceT>::masker_ptr();
    }
    return boost::python::extract<typename ImageSetData<BeamT, SequenceT>::masker_ptr>(masker)();
  }
  

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<MonochromaticBeam, Scan> > make_monochromatic_imageset_data1(boost::python::object reader,
                                                      boost::python::object masker) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<MonochromaticBeam, Scan> > self(
      new ImageSetData<MonochromaticBeam, Scan>(reader, make_masker_pointer<MonochromaticBeam, Scan>(masker)));

    // Return the imageset data
    return self;
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<MonochromaticBeam, Scan> > make_monochromatic_imageset_data2(boost::python::object reader,
                                                      boost::python::object masker,
                                                      std::string filename_template,
                                                      std::string vendor,
                                                      boost::python::dict params,
                                                      boost::python::object format) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<MonochromaticBeam, Scan> > self(
      new ImageSetData<MonochromaticBeam, Scan>(reader, make_masker_pointer<MonochromaticBeam, Scan>(masker)));

    // Set some stuff
    self->set_template(filename_template);
    self->set_vendor(vendor);
    self->set_params(detail::pickle_dumps(params));
    self->set_format(detail::pickle_dumps(format));

    // Return the imageset data
    return self;
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<TOFBeam, TOFSequence> > make_tof_imageset_data1(boost::python::object reader,
                                                      boost::python::object masker) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<TOFBeam, TOFSequence> > self(
      new ImageSetData<TOFBeam, TOFSequence>(reader, make_masker_pointer<TOFBeam, TOFSequence>(masker)));

    // Return the imageset data
    return self;
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<TOFBeam, TOFSequence> > make_tof_imageset_data2(boost::python::object reader,
                                                      boost::python::object masker,
                                                      std::string filename_template,
                                                      std::string vendor,
                                                      boost::python::dict params,
                                                      boost::python::object format) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<TOFBeam, TOFSequence> > self(
      new ImageSetData<TOFBeam, TOFSequence>(reader, make_masker_pointer<TOFBeam, TOFSequence>(masker)));

    // Set some stuff
    self->set_template(filename_template);
    self->set_vendor(vendor);
    self->set_params(detail::pickle_dumps(params));
    self->set_format(detail::pickle_dumps(format));

    // Return the imageset data
    return self;
  }

  /**
   * Get the parameters
   */
  template<typename BeamT, typename SequenceT>
  boost::python::object ImageSetData_get_params(ImageSetData<BeamT, SequenceT> &self) {
    return detail::pickle_loads(self.get_params());
  }

  /**
   * Set the parameters
   */
  template<typename BeamT, typename SequenceT>
  void ImageSetData_set_params(ImageSetData<BeamT, SequenceT> &self, boost::python::dict params) {
    self.set_params(detail::pickle_dumps(params));
  }

  /**
   * Get the format class
   */
  template<typename BeamT, typename SequenceT>
  boost::python::object ImageSetData_get_format(ImageSetData<BeamT, SequenceT> &self) {
    return detail::pickle_loads(self.get_format());
  }

  /**
   * Set the format class
   */

  template<typename BeamT, typename SequenceT>
  void ImageSetData_set_format(ImageSetData<BeamT, SequenceT> &self, boost::python::dict format) {
    self.set_format(detail::pickle_dumps(format));
  }

  /**
   * A constructor for the imageset class
   */
  typename boost::shared_ptr<ImageSet<MonochromaticBeam, Scan> > make_imageset(const ImageSetData<MonochromaticBeam, Scan> &data,
                                            boost::python::object indices) {
    if (indices == boost::python::object()) {
      return typename boost::shared_ptr<ImageSet<MonochromaticBeam, Scan> >(new ImageSet<MonochromaticBeam, Scan>(data));
    }

    return typename boost::shared_ptr<ImageSet<MonochromaticBeam, Scan> >(new ImageSet<MonochromaticBeam, Scan>(
      data, boost::python::extract<scitbx::af::const_ref<std::size_t> >(indices)()));
  }

  /**
   * Implement pickling for ImageSetData class
   */
  template<class BeamT, class SequenceT>
  struct ImageSetDataPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSetData<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(obj.reader(), obj.masker());
    }

    static boost::shared_ptr<BeamT> get_beam(const ImageSetData<BeamT, SequenceT> &self,
                                                std::size_t i) {
      return self.get_beam(i);
    }

    static boost::shared_ptr<Detector> get_detector(const ImageSetData<BeamT, SequenceT> &self,
                                                    std::size_t i) {
      return self.get_detector(i);
    }

    static boost::shared_ptr<Goniometer> get_goniometer(const ImageSetData<BeamT, SequenceT> &self,
                                                        std::size_t i) {
      return self.get_goniometer(i);
    }

    static boost::shared_ptr<SequenceT> get_sequence(const ImageSetData<BeamT, SequenceT> &self, std::size_t i) {
      return self.get_sequence(i);
    }

    template <typename Model, typename Func>
    static boost::python::tuple get_model_list(ImageSetData<BeamT, SequenceT> obj, Func get) {
      // Create a list of models and a list of indices
      std::vector<boost::shared_ptr<Model> > model_list;
      std::vector<std::size_t> index_list;
      for (std::size_t i = 0; i < obj.size(); ++i) {
        boost::shared_ptr<Model> m = get(obj, i);
        std::size_t k = model_list.size();
        for (std::size_t j = 0; j < k; ++j) {
          if (m.get() == model_list[j].get()) {
            k = j;
            break;
          }
        }
        if (k == model_list.size()) {
          model_list.push_back(m);
        }
        index_list.push_back(k);
      }

      // Convert to python lists and return tuple
      boost::python::list models;
      boost::python::list indices;
      for (std::size_t i = 0; i < model_list.size(); ++i) {
        models.append(model_list[i]);
      }
      for (std::size_t i = 0; i < index_list.size(); ++i) {
        indices.append(index_list[i]);
      }
      return boost::python::make_tuple(models, indices);
    }

    static boost::python::tuple get_model_tuple(ImageSetData<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(
        ImageSetDataPickleSuite::get_model_list<BeamT>(
          obj, &ImageSetDataPickleSuite::get_beam),
        ImageSetDataPickleSuite::get_model_list<Detector>(
          obj, &ImageSetDataPickleSuite::get_detector),
        ImageSetDataPickleSuite::get_model_list<Goniometer>(
          obj, &ImageSetDataPickleSuite::get_goniometer),
        ImageSetDataPickleSuite::get_model_list<SequenceT>(
          obj, &ImageSetDataPickleSuite::get_sequence));
    }

    static boost::python::tuple get_lookup_tuple(ImageSetData<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(
        boost::python::make_tuple(obj.external_lookup().mask().get_filename(),
                                  obj.external_lookup().mask().get_data()),
        boost::python::make_tuple(obj.external_lookup().gain().get_filename(),
                                  obj.external_lookup().gain().get_data()),
        boost::python::make_tuple(obj.external_lookup().pedestal().get_filename(),
                                  obj.external_lookup().pedestal().get_data()),
        boost::python::make_tuple(obj.external_lookup().dx().get_filename(),
                                  obj.external_lookup().dx().get_data()),
        boost::python::make_tuple(obj.external_lookup().dy().get_filename(),
                                  obj.external_lookup().dy().get_data()));
    }

    static boost::python::tuple getstate(ImageSetData<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(ImageSetDataPickleSuite::get_model_tuple(obj),
                                       ImageSetDataPickleSuite::get_lookup_tuple(obj),
                                       obj.get_template(),
                                       obj.get_vendor(),
                                       detail::bytes_from_std_string(obj.get_params()),
                                       detail::bytes_from_std_string(obj.get_format()));
    }

    template <typename Model, typename Func>
    static void set_model_list(ImageSetData<BeamT, SequenceT> &obj, boost::python::tuple data, Func set) {
      // Extract to python lists
      boost::python::list models =
        boost::python::extract<boost::python::list>(data[0])();
      boost::python::list indices =
        boost::python::extract<boost::python::list>(data[1])();

      // Convert to c++ vectors
      std::vector<boost::shared_ptr<Model> > model_list;
      std::vector<std::size_t> index_list;
      for (std::size_t i = 0; i < boost::python::len(models); ++i) {
        model_list.push_back(
          boost::python::extract<boost::shared_ptr<Model> >(models[i])());
      }
      for (std::size_t i = 0; i < boost::python::len(indices); ++i) {
        index_list.push_back(boost::python::extract<std::size_t>(indices[i])());
      }

      // Set the models
      DXTBX_ASSERT(index_list.size() == obj.size());
      for (std::size_t i = 0; i < index_list.size(); ++i) {
        DXTBX_ASSERT(index_list[i] < model_list.size());
        ((&obj)->*set)(model_list[index_list[i]], i);
      }
    }

    static void set_model_tuple(ImageSetData<BeamT, SequenceT> &obj, boost::python::tuple models) {
      DXTBX_ASSERT(boost::python::len(models) == 4);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_model_list<BeamT>(
        obj,
        boost::python::extract<boost::python::tuple>(models[0])(),
        &ImageSetData<BeamT, SequenceT>::set_beam);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_model_list<Detector>(
        obj,
        boost::python::extract<boost::python::tuple>(models[1]),
        &ImageSetData<BeamT, SequenceT>::set_detector);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_model_list<Goniometer>(
        obj,
        boost::python::extract<boost::python::tuple>(models[2]),
        &ImageSetData<BeamT, SequenceT>::set_goniometer);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_model_list<SequenceT>(
        obj,
        boost::python::extract<boost::python::tuple>(models[3]),
        &ImageSetData<BeamT, SequenceT>::set_sequence);
    }

    template <typename Data, typename Func>
    static void set_lookup_item(ImageSetData<BeamT, SequenceT> &obj,
                                boost::python::tuple lookup,
                                Func item) {
      DXTBX_ASSERT(boost::python::len(lookup) == 2);

      // Extract
      std::string filename = boost::python::extract<std::string>(lookup[0])();
      Data data = boost::python::extract<Data>(lookup[1])();

      // Set the filename and data
      ((&obj.external_lookup())->*item)().set_filename(filename);
      ((&obj.external_lookup())->*item)().set_data(data);
    }

    static void set_lookup_tuple(ImageSetData<BeamT, SequenceT> &obj, boost::python::tuple lookup) {
      DXTBX_ASSERT(boost::python::len(lookup) == 5);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_item<Image<bool> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[0])(),
        &ExternalLookup::mask);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[1])(),
        &ExternalLookup::gain);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[2])(),
        &ExternalLookup::pedestal);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[3])(),
        &ExternalLookup::dx);
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[4])(),
        &ExternalLookup::dy);
    }

    static void setstate(ImageSetData<BeamT, SequenceT> &obj, boost::python::tuple state) {
      DXTBX_ASSERT(boost::python::len(state) == 6);

      // Set the models
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_model_tuple(
        obj, boost::python::extract<boost::python::tuple>(state[0])());

      // Set the lookup
      ImageSetDataPickleSuite<BeamT, SequenceT>::set_lookup_tuple(
        obj, boost::python::extract<boost::python::tuple>(state[1])());

      // Set the properties
      obj.set_template(boost::python::extract<std::string>(state[2])());
      obj.set_vendor(boost::python::extract<std::string>(state[3])());
      obj.set_params(boost::python::extract<std::string>(state[4])());
      obj.set_format(boost::python::extract<std::string>(state[5])());
    }
  };

  /**
   * Implement pickling for ImageSet class
   */
  template<class BeamT, class SequenceT>
  struct ImageSetPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSet<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(obj.data(), obj.indices());
    }
  };

  /**
   * Implement pickling for ImageGrid class
   */
  struct ImageGridPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageGrid obj) {
      return boost::python::make_tuple(obj.data(), obj.indices(), obj.get_grid_size());
    }
  };

  /**
   * Implement pickling for ImageSequence class
   */
  template<class BeamT, class SequenceT>
  struct ImageSequencePickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSequence<BeamT, SequenceT> obj) {
      return boost::python::make_tuple(obj.data(),
                                       obj.indices(),
                                       obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_sequence());
    }
  };

  /**
   * Implement pickling for TOFImageSequence class
   */
  struct TOFImageSequencePickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(TOFImageSequence obj) {
      return boost::python::make_tuple(obj.data(),
                                       obj.indices(),
                                       obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_sequence());
    }
  };

  /**
   * Implement pickling for RotImageSequence class
   */
  struct RotImageSequencePickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(RotImageSequence obj) {
      return boost::python::make_tuple(obj.data(),
                                       obj.indices(),
                                       obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_sequence());
    }
  };

  /**
   * Get the external lookup data as a tuple of images
   */
  template <typename T>
  boost::python::object ExternalLookupItem_get_data(const ExternalLookupItem<T> &obj) {
    // Get the image data
    Image<T> data = obj.get_data();

    // If empty then return None
    if (data.empty()) {
      return boost::python::object();
    }

    // Otherwise, put image data in a list
    boost::python::list result;
    for (std::size_t i = 0; i < data.n_tiles(); ++i) {
      result.append(data.tile(i).data());
    }

    // Return the image - don't copy
#if __cplusplus > 199711L
    return std::move(result);
#else
    return result;
#endif
  }

  /**
   * Set the external lookup data
   */
  template <typename T>
  void ExternalLookupItem_set_data(ExternalLookupItem<T> &obj,
                                   boost::python::object item) {
    typedef typename scitbx::af::flex<T>::type flex_type;

    Image<T> data;

    // If item is not None
    if (item != boost::python::object()) {
      std::string name =
        boost::python::extract<std::string>(item.attr("__class__").attr("__name__"))();

      if (name == "tuple") {
        // If we have a tuple then add items of tuple to image data
        for (std::size_t i = 0; i < boost::python::len(item); ++i) {
          flex_type a = boost::python::extract<flex_type>(item)();
          data.push_back(ImageTile<T>(scitbx::af::versa<T, scitbx::af::c_grid<2> >(
            a.handle(), scitbx::af::c_grid<2>(a.accessor()))));
        }
      } else {
        try {
          // If we have a single array then add
          flex_type a = boost::python::extract<flex_type>(item)();
          data.push_back(ImageTile<T>(scitbx::af::versa<T, scitbx::af::c_grid<2> >(
            a.handle(), scitbx::af::c_grid<2>(a.accessor()))));

        } catch (boost::python::error_already_set) {
          data = boost::python::extract<Image<T> >(item)();
          boost::python::handle_exception();
        }
      }
    }

    // Set the image data
    obj.set_data(data);
  }

  template <typename T>
  boost::python::tuple image_as_tuple(const Image<T> &image) {
    boost::python::list result;
    for (std::size_t i = 0; i < image.n_tiles(); ++i) {
      result.append(image.tile(i).data());
    }
    return boost::python::tuple(result);
  }

  template<typename BeamT, typename SequenceT>
  boost::python::tuple ImageSet_get_raw_data(ImageSet<BeamT, SequenceT> &self, std::size_t index) {
    boost::python::tuple result;
    ImageBuffer buffer = self.get_raw_data(index);
    if (buffer.is_int()) {
      result = image_as_tuple<int>(buffer.as_int());
    } else if (buffer.is_double()) {
      result = image_as_tuple<double>(buffer.as_double());
    } else if (buffer.is_float()) {
      result = image_as_tuple<float>(buffer.as_float());
    } else {
      throw DXTBX_ERROR("Problem reading raw data");
    }
    return result;
  }

  template<typename BeamT, typename SequenceT>
  boost::python::tuple ImageSet_get_corrected_data(ImageSet<BeamT, SequenceT> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_corrected_data(index));
  }

  template<typename BeamT, typename SequenceT>
  boost::python::tuple ImageSet_get_gain(ImageSet<BeamT, SequenceT> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_gain(index));
  }

  template<typename BeamT, typename SequenceT>
  boost::python::tuple ImageSet_get_pedestal(ImageSet<BeamT, SequenceT> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_pedestal(index));
  }

  template<typename BeamT, typename SequenceT>
  boost::python::tuple ImageSet_get_mask(ImageSet<BeamT, SequenceT> &self, std::size_t index) {
    return image_as_tuple<bool>(self.get_mask(index));
  }

  /**
   * Wrapper for the external lookup items
   */
  template <typename T>
  void external_lookup_item_wrapper(const char *name) {
    using namespace boost::python;

    class_<ExternalLookupItem<T> >(name)
      .add_property("filename",
                    &ExternalLookupItem<T>::get_filename,
                    &ExternalLookupItem<T>::set_filename)
      .add_property(
        "data", &ExternalLookupItem<T>::get_data, &ExternalLookupItem<T>::set_data);
  }

  /**
   * If we have offset arrays set in the imageset then update the pixel to
   * millimeter strategy to use them
   */
  template<typename BeamT, typename SequenceT>
  void ImageSet_update_detector_px_mm_data(ImageSet<BeamT, SequenceT> &self) {
    Image<double> dx = self.external_lookup().dx().get_data();
    Image<double> dy = self.external_lookup().dy().get_data();
    DXTBX_ASSERT(dx.empty() == dy.empty());
    if (dx.empty() && dy.empty()) {
      return;
    }
    for (std::size_t i = 0; i < self.size(); ++i) {
      typename ImageSet<BeamT, SequenceT>::detector_ptr detector = self.get_detector_for_image(i);
      DXTBX_ASSERT(dx.n_tiles() == detector->size());
      DXTBX_ASSERT(dy.n_tiles() == detector->size());
      for (std::size_t i = 0; i < detector->size(); ++i) {
        Panel &panel = detector->operator[](i);
        if (panel.get_px_mm_strategy()->name() == "ParallaxCorrectedPxMmStrategy"
            || panel.get_px_mm_strategy()->name()
                 == "OffsetParallaxCorrectedPxMmStrategy") {
          boost::shared_ptr<OffsetParallaxCorrectedPxMmStrategy> strategy =
            boost::make_shared<OffsetParallaxCorrectedPxMmStrategy>(
              panel.get_mu(),
              panel.get_thickness(),
              dx.tile(i).data(),
              dy.tile(i).data());
          panel.set_px_mm_strategy(strategy);
        } else if (panel.get_px_mm_strategy()->name() == "SimplePxMmStrategy"
                   || panel.get_px_mm_strategy()->name() == "OffsetPxMmStrategy") {
          boost::shared_ptr<OffsetPxMmStrategy> strategy =
            boost::make_shared<OffsetPxMmStrategy>(dx.tile(i).data(),
                                                   dy.tile(i).data());
          panel.set_px_mm_strategy(strategy);
        }
      }
    }
  }

  /**
   * If we have offset arrays set in the imageset then update the pixel to
   * millimeter strategy to use them
   */
  template<typename BeamT, typename SequenceT>
  void ImageSequence_update_detector_px_mm_data(ImageSequence<BeamT, SequenceT> &self) {
    typename ImageSequence<BeamT, SequenceT>::detector_ptr detector = self.get_detector();
    Image<double> dx = self.external_lookup().dx().get_data();
    Image<double> dy = self.external_lookup().dy().get_data();
    DXTBX_ASSERT(dx.empty() == dy.empty());
    if (dx.empty() && dy.empty()) {
      return;
    }
    DXTBX_ASSERT(dx.n_tiles() == detector->size());
    DXTBX_ASSERT(dy.n_tiles() == detector->size());
    for (std::size_t i = 0; i < detector->size(); ++i) {
      Panel &panel = detector->operator[](i);
      if (panel.get_px_mm_strategy()->name() == "ParallaxCorrectedPxMmStrategy"
          || panel.get_px_mm_strategy()->name()
               == "OffsetParallaxCorrectedPxMmStrategy") {
        boost::shared_ptr<OffsetParallaxCorrectedPxMmStrategy> strategy =
          boost::make_shared<OffsetParallaxCorrectedPxMmStrategy>(panel.get_mu(),
                                                                  panel.get_thickness(),
                                                                  dx.tile(i).data(),
                                                                  dy.tile(i).data());
        panel.set_px_mm_strategy(strategy);
      } else if (panel.get_px_mm_strategy()->name() == "SimplePxMmStrategy"
                 || panel.get_px_mm_strategy()->name() == "OffsetPxMmStrategy") {
        boost::shared_ptr<OffsetPxMmStrategy> strategy =
          boost::make_shared<OffsetPxMmStrategy>(dx.tile(i).data(), dy.tile(i).data());
        panel.set_px_mm_strategy(strategy);
      }
    }
  }

  /**
   * Export the imageset classes
   */
  void export_imageset() {
    using namespace boost::python;

    external_lookup_item_wrapper<double>("ExternalLookupItemDouble");
    external_lookup_item_wrapper<bool>("ExternalLookupItemBool");

    class_<ExternalLookup>("ExternalLookup")
      .add_property("mask",
                    make_function(&ExternalLookup::mask, return_internal_reference<>()))
      .add_property("gain",
                    make_function(&ExternalLookup::gain, return_internal_reference<>()))
      .add_property(
        "pedestal",
        make_function(&ExternalLookup::pedestal, return_internal_reference<>()))
      .add_property("dx",
                    make_function(&ExternalLookup::dx, return_internal_reference<>()))
      .add_property("dy",
                    make_function(&ExternalLookup::dy, return_internal_reference<>()));

    class_<ImageSetData<TOFBeam, TOFSequence>, boost::shared_ptr<ImageSetData<TOFBeam, TOFSequence> > >("TOFImageSetData", no_init)
      .def("__init__",
           make_constructor(&make_tof_imageset_data1,
                            default_call_policies(),
                            (arg("reader"), arg("masker"))))
      .def("__init__",
           make_constructor(&make_tof_imageset_data2,
                            default_call_policies(),
                            (arg("reader"),
                             arg("masker"),
                             arg("template") = "",
                             arg("vendor") = "",
                             arg("params") = boost::python::object(),
                             arg("format") = boost::python::object())))
      .def("reader", &ImageSetData<TOFBeam, TOFSequence>::reader)
      .def("masker", &ImageSetData<TOFBeam, TOFSequence>::masker)
      .def("get_data", static_cast<ImageBuffer (ImageSetData<TOFBeam, TOFSequence>::*)(std::size_t)>(&ImageSetData<TOFBeam, TOFSequence>::get_data))
      .def("has_single_file_reader", &ImageSetData<TOFBeam, TOFSequence>::has_single_file_reader)
      .def("get_path", &ImageSetData<TOFBeam, TOFSequence>::get_path)
      .def("get_master_path", &ImageSetData<TOFBeam, TOFSequence>::get_master_path)
      .def("get_image_identifier", &ImageSetData<TOFBeam, TOFSequence>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetData<TOFBeam, TOFSequence>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetData<TOFBeam, TOFSequence>::is_marked_for_rejection)
      .def("get_beam", &ImageSetData<TOFBeam, TOFSequence>::get_beam)
      .def("get_detector", &ImageSetData<TOFBeam, TOFSequence>::get_detector)
      .def("get_goniometer", &ImageSetData<TOFBeam, TOFSequence>::get_goniometer)
      .def("get_sequence", &ImageSetData<TOFBeam, TOFSequence>::get_sequence)
      .def("set_beam", &ImageSetData<TOFBeam, TOFSequence>::set_beam)
      .def("set_detector", &ImageSetData<TOFBeam, TOFSequence>::set_detector)
      .def("set_goniometer", &ImageSetData<TOFBeam, TOFSequence>::set_goniometer)
      .def("set_sequence", &ImageSetData<TOFBeam, TOFSequence>::set_sequence)
      .def("get_template", &ImageSetData<TOFBeam, TOFSequence>::get_template)
      .def("set_template", &ImageSetData<TOFBeam, TOFSequence>::set_template)
      .def("get_vendor", &ImageSetData<TOFBeam, TOFSequence>::get_vendor)
      .def("set_vendor", &ImageSetData<TOFBeam, TOFSequence>::set_vendor)
      .def("get_params", &ImageSetData_get_params<TOFBeam, TOFSequence>)
      .def("set_params", &ImageSetData_set_params<TOFBeam, TOFSequence>)
      .def("get_format_class", &ImageSetData_get_format<TOFBeam, TOFSequence>)
      .def("set_format_class", &ImageSetData_set_format<TOFBeam, TOFSequence>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetData<TOFBeam, TOFSequence>::external_lookup, return_internal_reference<>()))
      .def_pickle(ImageSetDataPickleSuite<TOFBeam, TOFSequence>());

    class_<ImageSetData<MonochromaticBeam, Scan>, boost::shared_ptr<ImageSetData<MonochromaticBeam, Scan> > >("ImageSetData", no_init)
      .def("__init__",
           make_constructor(&make_monochromatic_imageset_data1,
                            default_call_policies(),
                            (arg("reader"), arg("masker"))))
      .def("__init__",
           make_constructor(&make_monochromatic_imageset_data2,
                            default_call_policies(),
                            (arg("reader"),
                             arg("masker"),
                             arg("template") = "",
                             arg("vendor") = "",
                             arg("params") = boost::python::object(),
                             arg("format") = boost::python::object())))
      .def("reader", &ImageSetData<MonochromaticBeam, Scan>::reader)
      .def("masker", &ImageSetData<MonochromaticBeam, Scan>::masker)
      .def("get_data", static_cast<ImageBuffer (ImageSetData<MonochromaticBeam, Scan>::*)(std::size_t)>(&ImageSetData<MonochromaticBeam, Scan>::get_data))
      .def("has_single_file_reader", &ImageSetData<MonochromaticBeam, Scan>::has_single_file_reader)
      .def("get_path", &ImageSetData<MonochromaticBeam, Scan>::get_path)
      .def("get_master_path", &ImageSetData<MonochromaticBeam, Scan>::get_master_path)
      .def("get_image_identifier", &ImageSetData<MonochromaticBeam, Scan>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetData<MonochromaticBeam, Scan>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetData<MonochromaticBeam, Scan>::is_marked_for_rejection)
      .def("get_beam", &ImageSetData<MonochromaticBeam, Scan>::get_beam)
      .def("get_detector", &ImageSetData<MonochromaticBeam, Scan>::get_detector)
      .def("get_goniometer", &ImageSetData<MonochromaticBeam, Scan>::get_goniometer)
      .def("get_sequence", &ImageSetData<MonochromaticBeam, Scan>::get_sequence)
      .def("set_beam", &ImageSetData<MonochromaticBeam, Scan>::set_beam)
      .def("set_detector", &ImageSetData<MonochromaticBeam, Scan>::set_detector)
      .def("set_goniometer", &ImageSetData<MonochromaticBeam, Scan>::set_goniometer)
      .def("set_sequence", &ImageSetData<MonochromaticBeam, Scan>::set_sequence)
      .def("get_template", &ImageSetData<MonochromaticBeam, Scan>::get_template)
      .def("set_template", &ImageSetData<MonochromaticBeam, Scan>::set_template)
      .def("get_vendor", &ImageSetData<MonochromaticBeam, Scan>::get_vendor)
      .def("set_vendor", &ImageSetData<MonochromaticBeam, Scan>::set_vendor)
      .def("get_params", &ImageSetData_get_params<MonochromaticBeam, Scan>)
      .def("set_params", &ImageSetData_set_params<MonochromaticBeam, Scan>)
      .def("get_format_class", &ImageSetData_get_format<MonochromaticBeam, Scan>)
      .def("set_format_class", &ImageSetData_set_format<MonochromaticBeam, Scan>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetData<MonochromaticBeam, Scan>::external_lookup, return_internal_reference<>()))
      .def_pickle(ImageSetDataPickleSuite<MonochromaticBeam, Scan>());

    class_<ImageSet<MonochromaticBeam, Scan> >("ImageSet", no_init)
      .def("__init__",
           make_constructor(&make_imageset,
                            default_call_policies(),
                            (arg("data"), arg("indices") = boost::python::object())))
      .def("data", &ImageSet<MonochromaticBeam, Scan>::data)
      .def("indices", &ImageSet<MonochromaticBeam, Scan>::indices)
      .def("size", &ImageSet<MonochromaticBeam, Scan>::size)
      .def("__len__", &ImageSet<MonochromaticBeam, Scan>::size)
      .def("has_dynamic_mask", &ImageSet<MonochromaticBeam, Scan>::has_dynamic_mask)
      .def("get_raw_data", &ImageSet_get_raw_data<MonochromaticBeam, Scan>)
      .def("get_corrected_data", &ImageSet_get_corrected_data<MonochromaticBeam, Scan>)
      .def("get_gain", &ImageSet_get_gain<MonochromaticBeam, Scan>)
      .def("get_pedestal", &ImageSet_get_pedestal<MonochromaticBeam, Scan>)
      .def("get_mask", &ImageSet_get_mask<MonochromaticBeam, Scan>)
      .def("get_beam", &ImageSet<MonochromaticBeam, Scan>::get_beam_for_image, (arg("index") = 0))
      .def("get_detector", &ImageSet<MonochromaticBeam, Scan>::get_detector_for_image, (arg("index") = 0))
      .def("get_goniometer", &ImageSet<MonochromaticBeam, Scan>::get_goniometer_for_image, (arg("index") = 0))
      .def("get_sequence", &ImageSet<MonochromaticBeam, Scan>::get_sequence_for_image, (arg("index") = 0))
      .def("set_beam", &ImageSet<MonochromaticBeam, Scan>::set_beam_for_image, (arg("index") = 0))
      .def("set_detector", &ImageSet<MonochromaticBeam, Scan>::set_detector_for_image, (arg("index") = 0))
      .def("set_goniometer", &ImageSet<MonochromaticBeam, Scan>::set_goniometer_for_image, (arg("index") = 0))
      .def("set_sequence", &ImageSet<MonochromaticBeam, Scan>::set_sequence_for_image, (arg("index") = 0))
      .def("get_path", &ImageSet<MonochromaticBeam, Scan>::get_path)
      .def("get_image_identifier", &ImageSet<MonochromaticBeam, Scan>::get_image_identifier)
      .def("mark_for_rejection", &ImageSet<MonochromaticBeam, Scan>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSet<MonochromaticBeam, Scan>::is_marked_for_rejection)
      .def("as_imageset", &ImageSet<MonochromaticBeam, Scan>::as_imageset)
      .def("complete_set", &ImageSet<MonochromaticBeam, Scan>::complete_set)
      .def("partial_set", &ImageSet<MonochromaticBeam, Scan>::partial_set)
      .def("clear_cache", &ImageSet<MonochromaticBeam, Scan>::clear_cache)
      .def("__eq__", &ImageSet<MonochromaticBeam, Scan>::operator==)
      .def("__ne__", &ImageSet<MonochromaticBeam, Scan>::operator!=)
      .def("update_detector_px_mm_data", &ImageSet_update_detector_px_mm_data<MonochromaticBeam, Scan>)
      .add_property(
        "external_lookup",
        make_function(&ImageSet<MonochromaticBeam, Scan>::external_lookup, return_internal_reference<>()))
      .def_pickle(ImageSetPickleSuite<MonochromaticBeam, Scan>());

    class_<ImageSet<TOFBeam, TOFSequence> >("TOFImageSet", no_init)
      .def("data", &ImageSet<TOFBeam, TOFSequence>::data)
      .def("indices", &ImageSet<TOFBeam, TOFSequence>::indices)
      .def("size", &ImageSet<TOFBeam, TOFSequence>::size)
      .def("__len__", &ImageSet<TOFBeam, TOFSequence>::size)
      .def("has_dynamic_mask", &ImageSet<TOFBeam, TOFSequence>::has_dynamic_mask)
      .def("get_raw_data", &ImageSet_get_raw_data<TOFBeam, TOFSequence>)
      .def("get_corrected_data", &ImageSet_get_corrected_data<TOFBeam, TOFSequence>)
      .def("get_gain", &ImageSet_get_gain<TOFBeam, TOFSequence>)
      .def("get_pedestal", &ImageSet_get_pedestal<TOFBeam, TOFSequence>)
      .def("get_mask", &ImageSet_get_mask<TOFBeam, TOFSequence>)
      .def("get_beam", &ImageSet<TOFBeam, TOFSequence>::get_beam_for_image, (arg("index") = 0))
      .def("get_detector", &ImageSet<TOFBeam, TOFSequence>::get_detector_for_image, (arg("index") = 0))
      .def("get_goniometer", &ImageSet<TOFBeam, TOFSequence>::get_goniometer_for_image, (arg("index") = 0))
      .def("get_sequence", &ImageSet<TOFBeam, TOFSequence>::get_sequence_for_image, (arg("index") = 0))
      .def("set_beam", &ImageSet<TOFBeam, TOFSequence>::set_beam_for_image, (arg("index") = 0))
      .def("set_detector", &ImageSet<TOFBeam, TOFSequence>::set_detector_for_image, (arg("index") = 0))
      .def("set_goniometer", &ImageSet<TOFBeam, TOFSequence>::set_goniometer_for_image, (arg("index") = 0))
      .def("set_sequence", &ImageSet<TOFBeam, TOFSequence>::set_sequence_for_image, (arg("index") = 0))
      .def("get_path", &ImageSet<TOFBeam, TOFSequence>::get_path)
      .def("get_image_identifier", &ImageSet<TOFBeam, TOFSequence>::get_image_identifier)
      .def("mark_for_rejection", &ImageSet<TOFBeam, TOFSequence>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSet<TOFBeam, TOFSequence>::is_marked_for_rejection)
      .def("as_imageset", &ImageSet<TOFBeam, TOFSequence>::as_imageset)
      .def("complete_set", &ImageSet<TOFBeam, TOFSequence>::complete_set)
      .def("partial_set", &ImageSet<TOFBeam, TOFSequence>::partial_set)
      .def("clear_cache", &ImageSet<TOFBeam, TOFSequence>::clear_cache)
      .def("__eq__", &ImageSet<TOFBeam, TOFSequence>::operator==)
      .def("__ne__", &ImageSet<TOFBeam, TOFSequence>::operator!=)
      .def("update_detector_px_mm_data", &ImageSet_update_detector_px_mm_data<TOFBeam, TOFSequence>)
      .add_property(
        "external_lookup",
        make_function(&ImageSet<TOFBeam, TOFSequence>::external_lookup, return_internal_reference<>()));

    class_<ImageGrid, bases<ImageSet<MonochromaticBeam, Scan> > >("ImageGrid", no_init)
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &, int2>((arg("data"), arg("grid_size"))))
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &, const scitbx::af::const_ref<std::size_t> &, int2>(
        (arg("data"), arg("indices"), arg("grid_size"))))
      .def("get_grid_size", &ImageGrid::get_grid_size)
      .def("from_imageset", &ImageGrid::from_imageset)
      .staticmethod("from_imageset")
      .def_pickle(ImageGridPickleSuite());

    class_<ImageSequence<MonochromaticBeam, Scan>, bases<ImageSet<MonochromaticBeam, Scan> > >("RotImageSequenceBase", no_init)
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &,
                const ImageSequence<MonochromaticBeam, Scan>::beam_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::detector_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::goniometer_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::sequence_ptr &>(
        (arg("data"), arg("beam"), arg("detector"), arg("goniometer"), arg("sequence"))))
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &,
                const scitbx::af::const_ref<std::size_t> &,
                const ImageSequence<MonochromaticBeam, Scan>::beam_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::detector_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::goniometer_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::sequence_ptr &>((arg("data"),
                                                  arg("indices"),
                                                  arg("beam"),
                                                  arg("detector"),
                                                  arg("goniometer"),
                                                  arg("sequence"))))
      .def("get_beam", &ImageSequence<MonochromaticBeam, Scan>::get_beam_for_image)
      .def("get_detector", &ImageSequence<MonochromaticBeam, Scan>::get_detector_for_image)
      .def("get_goniometer", &ImageSequence<MonochromaticBeam, Scan>::get_goniometer_for_image)
      .def("get_sequence", &ImageSequence<MonochromaticBeam, Scan>::get_sequence_for_image)
      .def("set_beam", &ImageSequence<MonochromaticBeam, Scan>::set_beam_for_image)
      .def("set_detector", &ImageSequence<MonochromaticBeam, Scan>::set_detector_for_image)
      .def("set_goniometer", &ImageSequence<MonochromaticBeam, Scan>::set_goniometer_for_image)
      .def("set_sequence", &ImageSequence<MonochromaticBeam, Scan>::set_sequence_for_image)
      .def("get_beam", &ImageSequence<MonochromaticBeam, Scan>::get_beam)
      .def("get_detector", &ImageSequence<MonochromaticBeam, Scan>::get_detector)
      .def("get_goniometer", &ImageSequence<MonochromaticBeam, Scan>::get_goniometer)
      .def("get_sequence", &ImageSequence<MonochromaticBeam, Scan>::get_sequence)
      .def("set_beam", &ImageSequence<MonochromaticBeam, Scan>::set_beam)
      .def("set_detector", &ImageSequence<MonochromaticBeam, Scan>::set_detector)
      .def("set_goniometer", &ImageSequence<MonochromaticBeam, Scan>::set_goniometer)
      .def("set_sequence", &ImageSequence<MonochromaticBeam, Scan>::set_sequence)
      .def("get_array_range", &ImageSequence<MonochromaticBeam, Scan>::get_array_range)
      .def("complete_set", &ImageSequence<MonochromaticBeam, Scan>::complete_sequence)
      .def("partial_set", &ImageSequence<MonochromaticBeam, Scan>::partial_sequence)
      .def("update_detector_px_mm_data", &ImageSequence_update_detector_px_mm_data<MonochromaticBeam, Scan>);


    class_<RotImageSequence, bases<ImageSequence<MonochromaticBeam, Scan> > >("RotImageSequence", no_init)
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &,
                const ImageSequence<MonochromaticBeam, Scan>::beam_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::detector_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::goniometer_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::sequence_ptr &>(
        (arg("data"), arg("beam"), arg("detector"), arg("goniometer"), arg("sequence"))))
      .def(init<const ImageSetData<MonochromaticBeam, Scan> &,
                const scitbx::af::const_ref<std::size_t> &,
                const ImageSequence<MonochromaticBeam, Scan>::beam_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::detector_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::goniometer_ptr &,
                const ImageSequence<MonochromaticBeam, Scan>::sequence_ptr &>((arg("data"),
                                                  arg("indices"),
                                                  arg("beam"),
                                                  arg("detector"),
                                                  arg("goniometer"),
                                                  arg("sequence"))))
      .def("partial_set", &RotImageSequence::partial_sequence)
      .def_pickle(RotImageSequencePickleSuite());

    class_<ImageSequence<TOFBeam, TOFSequence>, bases<ImageSet<TOFBeam, TOFSequence> > >("TOFImageSequenceBase", no_init)
      .def(init<const ImageSetData<TOFBeam, TOFSequence> &,
                const ImageSequence<TOFBeam, TOFSequence>::beam_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::detector_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::goniometer_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::sequence_ptr &>(
        (arg("data"), arg("beam"), arg("detector"), arg("goniometer"), arg("sequence"))))
      .def(init<const ImageSetData<TOFBeam, TOFSequence> &,
                const scitbx::af::const_ref<std::size_t> &,
                const ImageSequence<TOFBeam, TOFSequence>::beam_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::detector_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::goniometer_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::sequence_ptr &>((arg("data"),
                                                  arg("indices"),
                                                  arg("beam"),
                                                  arg("detector"),
                                                  arg("goniometer"),
                                                  arg("sequence"))))
      .def("get_beam", &ImageSequence<TOFBeam, TOFSequence>::get_beam_for_image)
      .def("get_detector", &ImageSequence<TOFBeam, TOFSequence>::get_detector_for_image)
      .def("get_goniometer", &ImageSequence<TOFBeam, TOFSequence>::get_goniometer_for_image)
      .def("get_sequence", &ImageSequence<TOFBeam, TOFSequence>::get_sequence_for_image)
      .def("set_beam", &ImageSequence<TOFBeam, TOFSequence>::set_beam_for_image)
      .def("set_detector", &ImageSequence<TOFBeam, TOFSequence>::set_detector_for_image)
      .def("set_goniometer", &ImageSequence<TOFBeam, TOFSequence>::set_goniometer_for_image)
      .def("set_sequence", &ImageSequence<TOFBeam, TOFSequence>::set_sequence_for_image)
      .def("get_beam", &ImageSequence<TOFBeam, TOFSequence>::get_beam)
      .def("get_detector", &ImageSequence<TOFBeam, TOFSequence>::get_detector)
      .def("get_goniometer", &ImageSequence<TOFBeam, TOFSequence>::get_goniometer)
      .def("get_sequence", &ImageSequence<TOFBeam, TOFSequence>::get_sequence)
      .def("set_beam", &ImageSequence<TOFBeam, TOFSequence>::set_beam)
      .def("set_detector", &ImageSequence<TOFBeam, TOFSequence>::set_detector)
      .def("set_goniometer", &ImageSequence<TOFBeam, TOFSequence>::set_goniometer)
      .def("set_sequence", &ImageSequence<TOFBeam, TOFSequence>::set_sequence)
      .def("get_array_range", &ImageSequence<TOFBeam, TOFSequence>::get_array_range)
      .def("complete_set", &ImageSequence<TOFBeam, TOFSequence>::complete_sequence)
      .def("partial_set", &ImageSequence<TOFBeam, TOFSequence>::partial_sequence)
      .def("update_detector_px_mm_data", &ImageSequence_update_detector_px_mm_data<TOFBeam, TOFSequence>);


    class_<TOFImageSequence, bases<ImageSequence<TOFBeam, TOFSequence> > >("TOFImageSequence", no_init)
      .def(init<const ImageSetData<TOFBeam, TOFSequence> &,
                const ImageSequence<TOFBeam, TOFSequence>::beam_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::detector_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::goniometer_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::sequence_ptr &>(
        (arg("data"), arg("beam"), arg("detector"), arg("goniometer"), arg("sequence"))))
      .def(init<const ImageSetData<TOFBeam, TOFSequence> &,
                const scitbx::af::const_ref<std::size_t> &,
                const ImageSequence<TOFBeam, TOFSequence>::beam_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::detector_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::goniometer_ptr &,
                const ImageSequence<TOFBeam, TOFSequence>::sequence_ptr &>((arg("data"),
                                                  arg("indices"),
                                                  arg("beam"),
                                                  arg("detector"),
                                                  arg("goniometer"),
                                                  arg("sequence"))))
      .def("partial_set", &TOFImageSequence::partial_sequence)
      .def_pickle(TOFImageSequencePickleSuite());

  }

  BOOST_PYTHON_MODULE(dxtbx_imageset_ext) {
      export_imageset();
  }

}}  // namespace dxtbx::boost_python
