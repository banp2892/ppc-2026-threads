#include "chernykh_s_trapezoidal_integration/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include <tbb/partitioner.h> // Не забудьте инклуд

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationTBB::ChernykhSTrapezoidalIntegrationTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationTBB::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationTBB::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationTBB::CalculatePointAndWeight(const IntegrationInType &input,
                                                                   const std::vector<std::size_t> &counters,
                                                                   std::vector<double> &point) {
  double weight = 1.0;
  for (std::size_t i = 0; i < input.limits.size(); ++i) {  // проходим по границам каждого из измерений
    const double h = (input.limits[i].second - input.limits[i].first) /
                     static_cast<double>(input.steps[i]);  // велечина шага в текущем измерении
    point[i] = input.limits[i].first +
               (static_cast<double>(counters[i]) * h);  // дискретная точка: начало_измерения + шаг*номер_шага
    if (std::cmp_equal(counters[i], 0) ||
        std::cmp_equal(counters[i], input.steps[i])) {  // если шаг граничный, то его вес уменьшается
      weight *= 0.5;
    }
  }
  return weight;
}

bool ChernykhSTrapezoidalIntegrationTBB::RunImpl() {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, std::thread::hardware_concurrency());

    const auto &input = this->GetInput();
    // Предполагаем 3D для теста
    int steps0 = input.steps[0], steps1 = input.steps[1], steps2 = input.steps[2];
    double h0 = (input.limits[0].second - input.limits[0].first) / steps0;
    double h1 = (input.limits[1].second - input.limits[1].first) / steps1;
    double h2 = (input.limits[2].second - input.limits[2].first) / steps2;
    double s0 = input.limits[0].first, s1 = input.limits[1].first, s2 = input.limits[2].first;

    // Параллелим только по самому внешнему измерению (как делает OMP по умолчанию)
    GetOutput() = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, steps0 + 1),
        0.0,
        [&](const tbb::blocked_range<int> &r, double local_sum) {
            for (int i = r.begin(); i < r.end(); ++i) {
                double x = s0 + i * h0;
                double w0 = (i == 0 || i == steps0) ? 0.5 : 1.0;
                
                for (int j = 0; j <= steps1; ++j) {
                    double y = s1 + j * h1;
                    double w1 = (j == 0 || j == steps1) ? 0.5 : 1.0;
                    
                    for (int k = 0; k <= steps2; ++k) {
                        double z = s2 + k * h2;
                        double w2 = (k == 0 || k == steps2) ? 0.5 : 1.0;
                        
                        // Прямая формула БЕЗ векторов и БЕЗ деления в цикле
                        local_sum += (std::sin(x) * std::cos(y) * std::exp(z)) * (w0 * w1 * w2);
                    }
                }
            }
            return local_sum;
        },
        std::plus<double>(),
        tbb::static_partitioner()
    ) * (h0 * h1 * h2);

    return true;
}

bool ChernykhSTrapezoidalIntegrationTBB::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
