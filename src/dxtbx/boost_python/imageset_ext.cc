#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/shared_ptr.hpp>
#include "Python.h"

#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <vector>
#include <dxtbx/imageset.h>
#include <dxtbx/model/beam.h>
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

  template<typename Beam>
  typename ImageSetData<Beam>::masker_ptr make_masker_pointer(boost::python::object masker) {
    if (masker == boost::python::object()) {
      return typename ImageSetData<Beam>::masker_ptr();
    }
    return boost::python::extract<typename ImageSetData<Beam>::masker_ptr>(masker)();
  }
  

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<MonochromaticBeam> > make_monochromatic_imageset_data1(boost::python::object reader,
                                                      boost::python::object masker) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<MonochromaticBeam> > self(
      new ImageSetData<MonochromaticBeam>(reader, make_masker_pointer<MonochromaticBeam>(masker)));

    // Return the imageset data
    return self;
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<MonochromaticBeam> > make_monochromatic_imageset_data2(boost::python::object reader,
                                                      boost::python::object masker,
                                                      std::string filename_template,
                                                      std::string vendor,
                                                      boost::python::dict params,
                                                      boost::python::object format) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<MonochromaticBeam> > self(
      new ImageSetData<MonochromaticBeam>(reader, make_masker_pointer<MonochromaticBeam>(masker)));

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
  boost::shared_ptr<ImageSetData<TOFBeam> > make_tof_imageset_data1(boost::python::object reader,
                                                      boost::python::object masker) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<TOFBeam> > self(
      new ImageSetData<TOFBeam>(reader, make_masker_pointer<TOFBeam>(masker)));

    // Return the imageset data
    return self;
  }

  /**
   * A constructor for the imageset data class
   */
  boost::shared_ptr<ImageSetData<TOFBeam> > make_tof_imageset_data2(boost::python::object reader,
                                                      boost::python::object masker,
                                                      std::string filename_template,
                                                      std::string vendor,
                                                      boost::python::dict params,
                                                      boost::python::object format) {
    // Create the pointer
    boost::shared_ptr<ImageSetData<TOFBeam> > self(
      new ImageSetData<TOFBeam>(reader, make_masker_pointer<TOFBeam>(masker)));

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
  template<typename Beam>
  boost::python::object ImageSetData_get_params(ImageSetData<Beam> &self) {
    return detail::pickle_loads(self.get_params());
  }

  /**
   * Set the parameters
   */
  template<typename Beam>
  void ImageSetData_set_params(ImageSetData<Beam> &self, boost::python::dict params) {
    self.set_params(detail::pickle_dumps(params));
  }

  /**
   * Get the format class
   */
  template<typename Beam>
  boost::python::object ImageSetData_get_format(ImageSetData<Beam> &self) {
    return detail::pickle_loads(self.get_format());
  }

  /**
   * Set the format class
   */

  template<typename Beam>
  void ImageSetData_set_format(ImageSetData<Beam> &self, boost::python::dict format) {
    self.set_format(detail::pickle_dumps(format));
  }

  /**
   * A constructor for the imageset class
   */
  typename boost::shared_ptr<ImageSet> make_imageset(const ImageSetData<MonochromaticBeam> &data,
                                            boost::python::object indices) {
    if (indices == boost::python::object()) {
      return typename boost::shared_ptr<ImageSet>(new ImageSet(data));
    }

    return typename boost::shared_ptr<ImageSet>(new ImageSet(
      data, boost::python::extract<scitbx::af::const_ref<std::size_t> >(indices)()));
  }

  /**
   * A constructor for the TOF imageset class
   */
  typename boost::shared_ptr<TOFImageSet> make_tof_imageset(const ImageSetData<TOFBeam> &data,
                                            boost::python::object tof_in_seconds,
                                            boost::python::object tof_indices,
                                            boost::python::object indices) {

    if (tof_indices != boost::python::object() && indices != boost::python::object()) {
      return typename boost::shared_ptr<TOFImageSet>(new TOFImageSet(data,
        boost::python::extract<scitbx::af::shared<double> >(tof_in_seconds)(),
        boost::python::extract<scitbx::af::const_ref<std::size_t> >(tof_indices)(),
          boost::python::extract<scitbx::af::const_ref<std::size_t> >(indices)()));
    }
    else if (tof_indices != boost::python::object()){
      return typename boost::shared_ptr<TOFImageSet>(new TOFImageSet(data,
        boost::python::extract<scitbx::af::shared<double> >(tof_in_seconds)(),
        boost::python::extract<scitbx::af::const_ref<std::size_t> >(tof_indices)()));      
    }

    return typename boost::shared_ptr<TOFImageSet>(new TOFImageSet(data,
      boost::python::extract<scitbx::af::shared<double> >(tof_in_seconds)()));
  }

  /**
   * Implement pickling for ImageSetData class
   */
  template<class Beam>
  struct ImageSetDataPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSetData<Beam> obj) {
      return boost::python::make_tuple(obj.reader(), obj.masker());
    }

    static boost::shared_ptr<Beam> get_beam(const ImageSetData<Beam> &self,
                                                std::size_t i) {
      return self.get_beam(i);
    }

    static boost::shared_ptr<Detector> get_detector(const ImageSetData<Beam> &self,
                                                    std::size_t i) {
      return self.get_detector(i);
    }

    static boost::shared_ptr<Goniometer> get_goniometer(const ImageSetData<Beam> &self,
                                                        std::size_t i) {
      return self.get_goniometer(i);
    }

    static boost::shared_ptr<Scan> get_scan(const ImageSetData<Beam> &self, std::size_t i) {
      return self.get_scan(i);
    }

    template <typename Model, typename Func>
    static boost::python::tuple get_model_list(ImageSetData<Beam> obj, Func get) {
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

    static boost::python::tuple get_model_tuple(ImageSetData<Beam> obj) {
      return boost::python::make_tuple(
        ImageSetDataPickleSuite::get_model_list<Beam>(
          obj, &ImageSetDataPickleSuite::get_beam),
        ImageSetDataPickleSuite::get_model_list<Detector>(
          obj, &ImageSetDataPickleSuite::get_detector),
        ImageSetDataPickleSuite::get_model_list<Goniometer>(
          obj, &ImageSetDataPickleSuite::get_goniometer),
        ImageSetDataPickleSuite::get_model_list<Scan>(
          obj, &ImageSetDataPickleSuite::get_scan));
    }

    static boost::python::tuple get_lookup_tuple(ImageSetData<Beam> obj) {
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

    static boost::python::tuple getstate(ImageSetData<Beam> obj) {
      return boost::python::make_tuple(ImageSetDataPickleSuite::get_model_tuple(obj),
                                       ImageSetDataPickleSuite::get_lookup_tuple(obj),
                                       obj.get_template(),
                                       obj.get_vendor(),
                                       detail::bytes_from_std_string(obj.get_params()),
                                       detail::bytes_from_std_string(obj.get_format()));
    }

    template <typename Model, typename Func>
    static void set_model_list(ImageSetData<Beam> &obj, boost::python::tuple data, Func set) {
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

    static void set_model_tuple(ImageSetData<Beam> &obj, boost::python::tuple models) {
      DXTBX_ASSERT(boost::python::len(models) == 4);
      ImageSetDataPickleSuite<Beam>::set_model_list<Beam>(
        obj,
        boost::python::extract<boost::python::tuple>(models[0])(),
        &ImageSetData<Beam>::set_beam);
      ImageSetDataPickleSuite<Beam>::set_model_list<Detector>(
        obj,
        boost::python::extract<boost::python::tuple>(models[1]),
        &ImageSetData<Beam>::set_detector);
      ImageSetDataPickleSuite<Beam>::set_model_list<Goniometer>(
        obj,
        boost::python::extract<boost::python::tuple>(models[2]),
        &ImageSetData<Beam>::set_goniometer);
      ImageSetDataPickleSuite<Beam>::set_model_list<Scan>(
        obj,
        boost::python::extract<boost::python::tuple>(models[3]),
        &ImageSetData<Beam>::set_scan);
    }

    template <typename Data, typename Func>
    static void set_lookup_item(ImageSetData<Beam> &obj,
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

    static void set_lookup_tuple(ImageSetData<Beam> &obj, boost::python::tuple lookup) {
      DXTBX_ASSERT(boost::python::len(lookup) == 5);
      ImageSetDataPickleSuite<Beam>::set_lookup_item<Image<bool> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[0])(),
        &ExternalLookup::mask);
      ImageSetDataPickleSuite<Beam>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[1])(),
        &ExternalLookup::gain);
      ImageSetDataPickleSuite<Beam>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[2])(),
        &ExternalLookup::pedestal);
      ImageSetDataPickleSuite<Beam>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[3])(),
        &ExternalLookup::dx);
      ImageSetDataPickleSuite<Beam>::set_lookup_item<Image<double> >(
        obj,
        boost::python::extract<boost::python::tuple>(lookup[4])(),
        &ExternalLookup::dy);
    }

    static void setstate(ImageSetData<Beam> &obj, boost::python::tuple state) {
      DXTBX_ASSERT(boost::python::len(state) == 6);

      // Set the models
      ImageSetDataPickleSuite<Beam>::set_model_tuple(
        obj, boost::python::extract<boost::python::tuple>(state[0])());

      // Set the lookup
      ImageSetDataPickleSuite<Beam>::set_lookup_tuple(
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
  struct ImageSetPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSet obj) {
      return boost::python::make_tuple(obj.data(), obj.indices());
    }
  };

  /**
   * Implement pickling for TOFImageSet class
   */
  
  struct TOFImageSetPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(TOFImageSet obj) {
      return boost::python::make_tuple(obj.data(), obj.tof_in_seconds(),
                                       obj.tof_indices(), obj.indices());
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
  struct ImageSequencePickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(ImageSequence obj) {
      return boost::python::make_tuple(obj.data(),
                                       obj.indices(),
                                       obj.get_beam(),
                                       obj.get_detector(),
                                       obj.get_goniometer(),
                                       obj.get_scan());
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

  template<typename Beam>
  boost::python::tuple ImageSetBase_get_raw_data(ImageSetBase<Beam> &self, std::size_t index) {
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

  template<typename Beam>
  boost::python::tuple ImageSetBase_get_corrected_data(ImageSetBase<Beam> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_corrected_data(index));
  }

  boost::python::tuple TOFImageSet_get_raw_data(TOFImageSet &self, std::size_t image_index, std::size_t tof_index) {
    boost::python::tuple result;
    ImageBuffer buffer = self.get_raw_data(image_index, tof_index);
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

  boost::python::tuple TOFImageSet_get_corrected_data(TOFImageSet &self, std::size_t image_index, std::size_t tof_index) {
    return image_as_tuple<double>(self.get_corrected_data(image_index, tof_index));
  }

  template<typename Beam>
  boost::python::tuple ImageSetBase_get_gain(ImageSetBase<Beam> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_gain(index));
  }

  template<typename Beam>
  boost::python::tuple ImageSetBase_get_pedestal(ImageSetBase<Beam> &self, std::size_t index) {
    return image_as_tuple<double>(self.get_pedestal(index));
  }

  template<typename Beam>
  boost::python::tuple ImageSetBase_get_mask(ImageSetBase<Beam> &self, std::size_t index) {
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
  template<typename Beam>
  void ImageSetBase_update_detector_px_mm_data(ImageSetBase<Beam> &self) {
    Image<double> dx = self.external_lookup().dx().get_data();
    Image<double> dy = self.external_lookup().dy().get_data();
    DXTBX_ASSERT(dx.empty() == dy.empty());
    if (dx.empty() && dy.empty()) {
      return;
    }
    for (std::size_t i = 0; i < self.size(); ++i) {
      typename ImageSetBase<Beam>::detector_ptr detector = self.get_detector_for_image(i);
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
  void ImageSequence_update_detector_px_mm_data(ImageSequence &self) {
    ImageSequence::detector_ptr detector = self.get_detector();
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

    class_<ImageSetData<TOFBeam>, boost::shared_ptr<ImageSetData<TOFBeam> > >("TOFImageSetData", no_init)
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
      .def("reader", &ImageSetData<TOFBeam>::reader)
      .def("masker", &ImageSetData<TOFBeam>::masker)
      .def("get_data", static_cast<ImageBuffer (ImageSetData<TOFBeam>::*)(std::size_t, std::size_t)>(&ImageSetData<TOFBeam>::get_data))
      .def("has_single_file_reader", &ImageSetData<TOFBeam>::has_single_file_reader)
      .def("get_path", &ImageSetData<TOFBeam>::get_path)
      .def("get_master_path", &ImageSetData<TOFBeam>::get_master_path)
      .def("get_image_identifier", &ImageSetData<TOFBeam>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetData<TOFBeam>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetData<TOFBeam>::is_marked_for_rejection)
      .def("get_beam", &ImageSetData<TOFBeam>::get_beam)
      .def("get_detector", &ImageSetData<TOFBeam>::get_detector)
      .def("get_goniometer", &ImageSetData<TOFBeam>::get_goniometer)
      .def("get_scan", &ImageSetData<TOFBeam>::get_scan)
      .def("set_beam", &ImageSetData<TOFBeam>::set_beam)
      .def("set_detector", &ImageSetData<TOFBeam>::set_detector)
      .def("set_goniometer", &ImageSetData<TOFBeam>::set_goniometer)
      .def("set_scan", &ImageSetData<TOFBeam>::set_scan)
      .def("get_template", &ImageSetData<TOFBeam>::get_template)
      .def("set_template", &ImageSetData<TOFBeam>::set_template)
      .def("get_vendor", &ImageSetData<TOFBeam>::get_vendor)
      .def("set_vendor", &ImageSetData<TOFBeam>::set_vendor)
      .def("get_params", &ImageSetData_get_params<TOFBeam>)
      .def("set_params", &ImageSetData_set_params<TOFBeam>)
      .def("get_format_class", &ImageSetData_get_format<TOFBeam>)
      .def("set_format_class", &ImageSetData_set_format<TOFBeam>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetData<TOFBeam>::external_lookup, return_internal_reference<>()))
      .def_pickle(ImageSetDataPickleSuite<TOFBeam>());

    //.def("get_data", &ImageSetData<MonochromaticBeam>::get_data, arg("index"))
    class_<ImageSetData<MonochromaticBeam>, boost::shared_ptr<ImageSetData<MonochromaticBeam> > >("ImageSetData", no_init)
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
      .def("reader", &ImageSetData<MonochromaticBeam>::reader)
      .def("masker", &ImageSetData<MonochromaticBeam>::masker)
      .def("get_data", static_cast<ImageBuffer (ImageSetData<MonochromaticBeam>::*)(std::size_t)>(&ImageSetData<MonochromaticBeam>::get_data))
      .def("has_single_file_reader", &ImageSetData<MonochromaticBeam>::has_single_file_reader)
      .def("get_path", &ImageSetData<MonochromaticBeam>::get_path)
      .def("get_master_path", &ImageSetData<MonochromaticBeam>::get_master_path)
      .def("get_image_identifier", &ImageSetData<MonochromaticBeam>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetData<MonochromaticBeam>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetData<MonochromaticBeam>::is_marked_for_rejection)
      .def("get_beam", &ImageSetData<MonochromaticBeam>::get_beam)
      .def("get_detector", &ImageSetData<MonochromaticBeam>::get_detector)
      .def("get_goniometer", &ImageSetData<MonochromaticBeam>::get_goniometer)
      .def("get_scan", &ImageSetData<MonochromaticBeam>::get_scan)
      .def("set_beam", &ImageSetData<MonochromaticBeam>::set_beam)
      .def("set_detector", &ImageSetData<MonochromaticBeam>::set_detector)
      .def("set_goniometer", &ImageSetData<MonochromaticBeam>::set_goniometer)
      .def("set_scan", &ImageSetData<MonochromaticBeam>::set_scan)
      .def("get_template", &ImageSetData<MonochromaticBeam>::get_template)
      .def("set_template", &ImageSetData<MonochromaticBeam>::set_template)
      .def("get_vendor", &ImageSetData<MonochromaticBeam>::get_vendor)
      .def("set_vendor", &ImageSetData<MonochromaticBeam>::set_vendor)
      .def("get_params", &ImageSetData_get_params<MonochromaticBeam>)
      .def("set_params", &ImageSetData_set_params<MonochromaticBeam>)
      .def("get_format_class", &ImageSetData_get_format<MonochromaticBeam>)
      .def("set_format_class", &ImageSetData_set_format<MonochromaticBeam>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetData<MonochromaticBeam>::external_lookup, return_internal_reference<>()))
      .def_pickle(ImageSetDataPickleSuite<MonochromaticBeam>());




    class_<ImageSetBase<MonochromaticBeam>>("ImageSetBase", no_init)
      .def("data", &ImageSetBase<MonochromaticBeam>::data)
      .def("indices", &ImageSetBase<MonochromaticBeam>::indices)
      .def("size", &ImageSetBase<MonochromaticBeam>::size)
      .def("__len__", &ImageSetBase<MonochromaticBeam>::size)
      .def("has_dynamic_mask", &ImageSetBase<MonochromaticBeam>::has_dynamic_mask)
      .def("get_raw_data", &ImageSetBase_get_raw_data<MonochromaticBeam>)
      .def("get_corrected_data", &ImageSetBase_get_corrected_data<MonochromaticBeam>)
      .def("get_gain", &ImageSetBase_get_gain<MonochromaticBeam>)
      .def("get_pedestal", &ImageSetBase_get_pedestal<MonochromaticBeam>)
      .def("get_mask", &ImageSetBase_get_mask<MonochromaticBeam>)
      .def("get_beam", &ImageSetBase<MonochromaticBeam>::get_beam_for_image, (arg("index") = 0))
      .def("get_detector", &ImageSetBase<MonochromaticBeam>::get_detector_for_image, (arg("index") = 0))
      .def("get_goniometer", &ImageSetBase<MonochromaticBeam>::get_goniometer_for_image, (arg("index") = 0))
      .def("get_scan", &ImageSetBase<MonochromaticBeam>::get_scan_for_image, (arg("index") = 0))
      .def("set_beam", &ImageSetBase<MonochromaticBeam>::set_beam_for_image, (arg("index") = 0))
      .def("set_detector", &ImageSetBase<MonochromaticBeam>::set_detector_for_image, (arg("index") = 0))
      .def("set_goniometer", &ImageSetBase<MonochromaticBeam>::set_goniometer_for_image, (arg("index") = 0))
      .def("set_scan", &ImageSetBase<MonochromaticBeam>::set_scan_for_image, (arg("index") = 0))
      .def("get_path", &ImageSetBase<MonochromaticBeam>::get_path)
      .def("get_image_identifier", &ImageSetBase<MonochromaticBeam>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetBase<MonochromaticBeam>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetBase<MonochromaticBeam>::is_marked_for_rejection)
      .def("as_imageset", &ImageSetBase<MonochromaticBeam>::as_imageset)
      .def("complete_set", &ImageSetBase<MonochromaticBeam>::complete_set)
      .def("partial_set", &ImageSetBase<MonochromaticBeam>::partial_set)
      .def("clear_cache", &ImageSetBase<MonochromaticBeam>::clear_cache)
      .def("__eq__", &ImageSetBase<MonochromaticBeam>::operator==)
      .def("__ne__", &ImageSetBase<MonochromaticBeam>::operator!=)
      .def("update_detector_px_mm_data", &ImageSetBase_update_detector_px_mm_data<MonochromaticBeam>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetBase<MonochromaticBeam>::external_lookup, return_internal_reference<>()));

    class_<ImageSetBase<TOFBeam>>("ImageSetBase", no_init)
      .def("data", &ImageSetBase<TOFBeam>::data)
      .def("indices", &ImageSetBase<TOFBeam>::indices)
      .def("size", &ImageSetBase<TOFBeam>::size)
      .def("__len__", &ImageSetBase<TOFBeam>::size)
      .def("has_dynamic_mask", &ImageSetBase<TOFBeam>::has_dynamic_mask)
      .def("get_raw_data", &ImageSetBase_get_raw_data<TOFBeam>)
      .def("get_corrected_data", &ImageSetBase_get_corrected_data<TOFBeam>)
      .def("get_gain", &ImageSetBase_get_gain<TOFBeam>)
      .def("get_pedestal", &ImageSetBase_get_pedestal<TOFBeam>)
      .def("get_mask", &ImageSetBase_get_mask<TOFBeam>)
      .def("get_beam", &ImageSetBase<TOFBeam>::get_beam_for_image, (arg("index") = 0))
      .def("get_detector", &ImageSetBase<TOFBeam>::get_detector_for_image, (arg("index") = 0))
      .def("get_goniometer", &ImageSetBase<TOFBeam>::get_goniometer_for_image, (arg("index") = 0))
      .def("get_scan", &ImageSetBase<TOFBeam>::get_scan_for_image, (arg("index") = 0))
      .def("set_beam", &ImageSetBase<TOFBeam>::set_beam_for_image, (arg("index") = 0))
      .def("set_detector", &ImageSetBase<TOFBeam>::set_detector_for_image, (arg("index") = 0))
      .def("set_goniometer", &ImageSetBase<TOFBeam>::set_goniometer_for_image, (arg("index") = 0))
      .def("set_scan", &ImageSetBase<TOFBeam>::set_scan_for_image, (arg("index") = 0))
      .def("get_path", &ImageSetBase<TOFBeam>::get_path)
      .def("get_image_identifier", &ImageSetBase<TOFBeam>::get_image_identifier)
      .def("mark_for_rejection", &ImageSetBase<TOFBeam>::mark_for_rejection)
      .def("is_marked_for_rejection", &ImageSetBase<TOFBeam>::is_marked_for_rejection)
      .def("as_imageset", &ImageSetBase<TOFBeam>::as_imageset)
      .def("complete_set", &ImageSetBase<TOFBeam>::complete_set)
      .def("partial_set", &ImageSetBase<TOFBeam>::partial_set)
      .def("clear_cache", &ImageSetBase<TOFBeam>::clear_cache)
      .def("__eq__", &ImageSetBase<TOFBeam>::operator==)
      .def("__ne__", &ImageSetBase<TOFBeam>::operator!=)
      .def("update_detector_px_mm_data", &ImageSetBase_update_detector_px_mm_data<TOFBeam>)
      .add_property(
        "external_lookup",
        make_function(&ImageSetBase<TOFBeam>::external_lookup, return_internal_reference<>()));

    class_<ImageSet, bases<ImageSetBase<MonochromaticBeam> > >("ImageSet", no_init)
      .def("__init__",
           make_constructor(&make_imageset,
                            default_call_policies(),
                            (arg("data"), arg("indices") = boost::python::object())))
      .def_pickle(ImageSetPickleSuite());

    class_<TOFImageSet, bases<ImageSetBase<TOFBeam> > >("TOFImageSet", no_init)
      .def("__init__",
           make_constructor(&make_tof_imageset,
                            default_call_policies(),
                            (arg("data"), arg("tof_in_seconds") = boost::python::object(),
                            arg("tof_indices") = boost::python::object(),
                            arg("indices") = boost::python::object())))
      .def("get_corrected_data", &TOFImageSet_get_corrected_data)
      .def("get_raw_data", &TOFImageSet_get_raw_data)
      .def("partial_set", &TOFImageSet::partial_set)
      .def("tof_in_seconds", &TOFImageSet::tof_in_seconds)
      .def_pickle(TOFImageSetPickleSuite());
    


    class_<ImageGrid, bases<ImageSetBase<MonochromaticBeam> > >("ImageGrid", no_init)
      .def(init<const ImageSetData<MonochromaticBeam> &, int2>((arg("data"), arg("grid_size"))))
      .def(init<const ImageSetData<MonochromaticBeam> &, const scitbx::af::const_ref<std::size_t> &, int2>(
        (arg("data"), arg("indices"), arg("grid_size"))))
      .def("get_grid_size", &ImageGrid::get_grid_size)
      .def("from_imageset", &ImageGrid::from_imageset)
      .staticmethod("from_imageset")
      .def_pickle(ImageGridPickleSuite());

    class_<ImageSequence, bases<ImageSetBase<MonochromaticBeam> > >("ImageSequence", no_init)
      .def(init<const ImageSetData<MonochromaticBeam> &,
                const ImageSequence::beam_ptr &,
                const ImageSequence::detector_ptr &,
                const ImageSequence::goniometer_ptr &,
                const ImageSequence::scan_ptr &>(
        (arg("data"), arg("beam"), arg("detector"), arg("goniometer"), arg("scan"))))
      .def(init<const ImageSetData<MonochromaticBeam> &,
                const scitbx::af::const_ref<std::size_t> &,
                const ImageSequence::beam_ptr &,
                const ImageSequence::detector_ptr &,
                const ImageSequence::goniometer_ptr &,
                const ImageSequence::scan_ptr &>((arg("data"),
                                                  arg("indices"),
                                                  arg("beam"),
                                                  arg("detector"),
                                                  arg("goniometer"),
                                                  arg("scan"))))
      .def("get_beam", &ImageSequence::get_beam_for_image)
      .def("get_detector", &ImageSequence::get_detector_for_image)
      .def("get_goniometer", &ImageSequence::get_goniometer_for_image)
      .def("get_scan", &ImageSequence::get_scan_for_image)
      .def("set_beam", &ImageSequence::set_beam_for_image)
      .def("set_detector", &ImageSequence::set_detector_for_image)
      .def("set_goniometer", &ImageSequence::set_goniometer_for_image)
      .def("set_scan", &ImageSequence::set_scan_for_image)
      .def("get_beam", &ImageSequence::get_beam)
      .def("get_detector", &ImageSequence::get_detector)
      .def("get_goniometer", &ImageSequence::get_goniometer)
      .def("get_scan", &ImageSequence::get_scan)
      .def("set_beam", &ImageSequence::set_beam)
      .def("set_detector", &ImageSequence::set_detector)
      .def("set_goniometer", &ImageSequence::set_goniometer)
      .def("set_scan", &ImageSequence::set_scan)
      .def("get_array_range", &ImageSequence::get_array_range)
      .def("complete_set", &ImageSequence::complete_sequence)
      .def("partial_set", &ImageSequence::partial_sequence)
      .def("update_detector_px_mm_data", &ImageSequence_update_detector_px_mm_data)
      .def_pickle(ImageSequencePickleSuite());
  }

  BOOST_PYTHON_MODULE(dxtbx_imageset_ext) {
    export_imageset();
  }

}}  // namespace dxtbx::boost_python
