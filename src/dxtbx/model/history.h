#ifndef DXTBX_MODEL_HISTORY_H
#define DXTBX_MODEL_HISTORY_H

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/extract.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace dxtbx { namespace model {

  /**
   * This class keeps track of the serialization history of an experiment.
   */
  class History {
  public:
    History() {}

    History(const boost::python::list &history) {
      set_history_from_list(history);
    }

    void set_history(const std::vector<std::string> &history) {
      history_ = history;
    }

    void append_history_item(const std::string &dispatcher,
                             const std::string &version,
                             const std::string &flag) {
      // Timestamp as current UTC time
      boost::posix_time::ptime now_utc =
        boost::posix_time::second_clock::universal_time();
      std::string utc_string = boost::posix_time::to_iso_extended_string(now_utc) + "Z";

      // Format the message
      std::string message = utc_string + "|" + dispatcher + "|" + version;
      if (!flag.empty()) {
        message += "|" + flag;
      }
      history_.push_back(message);
    }

    void set_history_from_list(const boost::python::list &history) {
      history_.clear();

      long length = boost::python::len(history);
      history_.reserve(length);

      for (long i = 0; i < length; ++i) {
        boost::python::extract<std::string> extractor(history[i]);
        history_.push_back(extractor());
      }
    }

    std::vector<std::string> get_history() const {
      return history_;
    }

    boost::python::list get_history_as_list() const {
      boost::python::list result;
      for (const auto &item : history_) {
        result.append(item);
      }
      return result;
    }

  private:
    std::vector<std::string> history_;
  };

}}  // namespace dxtbx::model

#endif  // DXTBX_MODEL_HISTORY_H