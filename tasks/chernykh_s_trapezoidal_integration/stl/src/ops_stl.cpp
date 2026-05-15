#include "chernykh_s_trapezoidal_integration/stl/include/ops_stl.hpp"


#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>
#include <thread>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"

namespace chernykh_s_trapezoidal_integration {

ChernykhSTrapezoidalIntegrationSTL::ChernykhSTrapezoidalIntegrationSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ChernykhSTrapezoidalIntegrationSTL::ValidationImpl() {
  const auto &input = this->GetInput();
  if (input.limits.empty() || input.limits.size() != input.steps.size()) {
    return false;
  }
  return std::ranges::all_of(input.steps, [](int s) { return s > 0; });
}

bool ChernykhSTrapezoidalIntegrationSTL::PreProcessingImpl() {
  return true;
}

double ChernykhSTrapezoidalIntegrationSTL::CalculatePointAndWeight(const IntegrationInType &input,
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





bool ChernykhSTrapezoidalIntegrationSTL::RunImpl() {
  auto & input= GetInput();
  int dims = input.limits.size(); // количество измерений
  int64_t total_point = 1;
  for(int setka: input.steps){
    total_point *= static_cast<int64_t>(setka)+1;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads ==0) num_threads=1;
  std::cout<<num_threads<<std::endl;

  std::vector<double> result(num_threads, 0.0); // вектор, в котором будут храниться результаты
  std::vector<std::thread> threads(num_threads); // вектор потоков 
  std::vector<std::pair<int64_t,int64_t>> borders(num_threads);

  int64_t points_per_thread = total_point / num_threads;
  int64_t remainder = total_point % num_threads;

  int64_t start = 0;
  for(int i=0; i <num_threads;i++){ // 11 3 по 3 остаток 2 : [0,3+1), [4,8), [8,11)
    borders[i].first = start;
    start+=points_per_thread;
    if(i<remainder){
      start++;
    }
    borders[i].second = start;
  }
  auto WorkFunction = [&](int64_t start, int64_t end, double& local_result){
    local_result = 0.0;
    std::vector<double> point(dims);
    for(int64_t i = start; i < end;i++){
      // Метод трапеций : вес(точка)*f(точка)
      int64_t local_index = i;
      std::vector<size_t> grid_index(dims,0); // для хранения идексов точек
      
      for( int j = dims-1; j>=0;j--){ // считаем индекс
        grid_index[j] = local_index % (input.steps[j]+1);
        local_index /= input.steps[j]+1;
      }
      double weight = CalculatePointAndWeight(input, grid_index, point);
      local_result += input.func(point) * weight;
    }
  };
  for (int i = 0; i < num_threads; ++i) {
    threads[i] = std::thread(WorkFunction, borders[i].first, borders[i].second, std::ref(result[i]));
  }

  
for (int i = 0; i < num_threads; ++i) {
    threads[i].join();
  }
  double output_result = 0.0;
  for(double loc_res: result){
output_result+=loc_res;
  }

  for(int i=0;i<dims;i++){
    double h = (input.limits[i].second - input.limits[i].first) / input.steps[i];
    output_result*=h;
    }
  
  GetOutput() = output_result;
  return true;
}

bool ChernykhSTrapezoidalIntegrationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chernykh_s_trapezoidal_integration
