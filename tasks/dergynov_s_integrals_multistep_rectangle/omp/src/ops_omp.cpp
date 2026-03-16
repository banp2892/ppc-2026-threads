#include "dergynov_s_integrals_multistep_rectangle/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  for (const auto &[left, right] : borders) {
    if (!std::isfinite(left) || !std::isfinite(right)) return false;
    if (left >= right) return false;
  }
  return true;
}

std::vector<int> LinearToMultiIndex(size_t linear_idx, int dim, int n) {
  std::vector<int> idx(dim);
  size_t tmp = linear_idx;
  for (int d = dim - 1; d >= 0; --d) {
    idx[d] = tmp % n;
    tmp /= n;
  }
  return idx;
}

}  // namespace

DergynovSIntegralsMultistepRectangleOMP::DergynovSIntegralsMultistepRectangleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleOMP::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();

  if (!func) return false;
  if (n <= 0) return false;
  if (borders.empty()) return false;

  return ValidateBorders(borders);
}

bool DergynovSIntegralsMultistepRectangleOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleOMP::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;

  for (int i = 0; i < dim; ++i) {
    const double left = borders[i].first;
    const double right = borders[i].second;
    h[i] = (right - left) / n;
    cell_volume *= h[i];
  }

  size_t total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= n;
  }

  std::vector<double> local_sums(omp_get_max_threads(), 0.0);
  
  int error_flag = 0;

#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    double local_sum = 0.0;

#pragma omp for schedule(static)
    for (size_t linear_idx = 0; linear_idx < total_points; ++linear_idx) {
      if (error_flag) continue;

      std::vector<int> idx = LinearToMultiIndex(linear_idx, dim, n);

      std::vector<double> point(dim);
      for (int d = 0; d < dim; ++d) {
        point[d] = borders[d].first + (idx[d] + 0.5) * h[d];
      }

      double f_val = func(point);
      if (!std::isfinite(f_val)) {
        #pragma omp atomic write
        error_flag = 1;
        continue;
      }
      local_sum += f_val;
    }

    local_sums[thread_id] = local_sum;
  }

  if (error_flag) {
    return false;
  }

  double total_sum = 0.0;
  for (double s : local_sums) {
    total_sum += s;
  }

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
