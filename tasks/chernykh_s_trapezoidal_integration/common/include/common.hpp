#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace chernykh_s_trapezoidal_integration {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chernykh_s_trapezoidal_integration
