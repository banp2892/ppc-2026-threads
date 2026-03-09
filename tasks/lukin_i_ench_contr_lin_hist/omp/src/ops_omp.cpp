#include "lukin_i_ench_contr_lin_hist/omp/include/ops_omp.hpp"

#include <algorithm>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"

namespace lukin_i_ench_contr_lin_hist {

LukinITestTaskOMP::LukinITestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskOMP::ValidationImpl() {
  return !(GetInput().empty());
}

bool LukinITestTaskOMP::PreProcessingImpl() {
  return true;
}

bool LukinITestTaskOMP::RunImpl() {
  auto min_it = std::ranges::min_element(GetInput().begin(), GetInput().end());
  auto max_it = std::ranges::max_element(GetInput().begin(), GetInput().end());

  unsigned char min = *min_it;
  unsigned char max = *max_it;

  if (max == min)  // Однотонное изображение
  {
    GetOutput() = GetInput();
    return true;
  }

  float scale = 255.0F / static_cast<float>(max - min);

  int size = static_cast<int>(GetInput().size());

  for (int i = 0; i < size; i++) {  // Линейное растяжение
    GetOutput()[i] = static_cast<unsigned char>(static_cast<float>(GetInput()[i] - min) * scale);
  }

  return true;
}

bool LukinITestTaskOMP::PostProcessingImpl() {
  return true;
}

}  // namespace lukin_i_ench_contr_lin_hist
