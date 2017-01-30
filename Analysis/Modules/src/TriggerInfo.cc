#include "Modules/interface/TriggerInfo.h"
#include "Utilities/interface/JsonTools.h"
#include "Objects/interface/Unhash.h"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "boost/lexical_cast.hpp"
#include <fstream>

namespace ic {

TriggerInfo::TriggerInfo(std::string const& name) : ModuleBase(name) {
  triggerpaths_label_ = "triggerPaths";
  eventinfo_label_ = "eventInfo";
  output_file_ = "trigger_info.json";
}

TriggerInfo::~TriggerInfo() { ; }

int TriggerInfo::PreAnalysis() {
  PrintHeader("TriggerInfo");
  PrintArg("eventinfo_label", eventinfo_label_);
  PrintArg("triggerpaths_label", triggerpaths_label_);
  PrintArg("output_file", output_file_);
  return 0;
}

int TriggerInfo::Execute(TreeEvent* event) {
  EventInfo const* evt_info = event->GetPtr<EventInfo>(eventinfo_label_);
  int run = evt_info->run();
  auto const& trig_info = event->GetPtrVec<TriggerPath>(triggerpaths_label_);

  for (auto trg : trig_info) {
    std::string name = trg->name() != "" ? trg->name() : Unhash::Get(trg->id());
    unsigned version = trg->version();
    auto & info = info_[name][version];
    if (info.first_run < 0 || run < info.first_run) {
      info.first_run = run;
    }
    if (info.last_run < 0 || run > info.last_run) {
      info.last_run = run;
    }

    if (!info.yield_per_run.count(run)) {
      info.yield_per_run[run] = 0;
    }
    if (trg->accept()) {
      ++(info.yield_per_run[run]);
    }
  }
  return 0;
}

int TriggerInfo::PostAnalysis() {
  for (auto const& path_info : info_) {
    for (auto const& version_info : path_info.second) {
      Json::Value& details =
          output_["info"][path_info.first]
                 [boost::lexical_cast<std::string>(version_info.first)];
      details["first_run"] = version_info.second.first_run;
      details["last_run"] = version_info.second.last_run;
      for (auto const& yield_info : version_info.second.yield_per_run) {
        output_["yields"][path_info.first][boost::lexical_cast<std::string>(yield_info.first)] = yield_info.second;
      }
    }
  }
  std::ofstream output;
  output.open(output_file_);
  Json::StyledWriter writer;
  output << writer.write(output_);
  output.close();
  return 0;
}

}