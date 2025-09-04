#include <Rcpp.h>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

struct Grid2D {
  std::vector<std::vector<std::vector<int>>> cells;
  const NumericVector& mz_values;
  const NumericVector& rt_values;
  double mz_min, mz_max, rt_min, rt_max;
  int mz_bins, rt_bins;

  Grid2D(int mz_b, double mz_mn, double mz_mx,
         int rt_b, double rt_mn, double rt_mx,
         const NumericVector& mz_vals,
         const NumericVector& rt_vals)
    : mz_values(mz_vals), rt_values(rt_vals),
      mz_min(mz_mn), mz_max(mz_mx), rt_min(rt_mn), rt_max(rt_mx),
      mz_bins(mz_b), rt_bins(rt_b) {
    cells.resize(mz_bins, std::vector<std::vector<int>>(rt_bins));
  }

  void add_point(int idx, double mz, double rt) {
    int mz_idx = std::min(
      static_cast<int>((mz - mz_min) / (mz_max - mz_min) * mz_bins),
      mz_bins - 1
    );
    int rt_idx = std::min(
      static_cast<int>((rt - rt_min) / (rt_max - rt_min) * rt_bins),
      rt_bins - 1
    );
    mz_idx = std::max(0, mz_idx);
    rt_idx = std::max(0, rt_idx);
    cells[mz_idx][rt_idx].push_back(idx);
  }

  std::vector<int> get_nearby_points(
      double mz, double rt,
      double mz_window_1x, double rt_window_1x,
      bool use_2x_window = false
  ) const {
    std::vector<int> result;

    double mz_window = use_2x_window ? 2 * mz_window_1x : mz_window_1x;
    double rt_window = use_2x_window ? 2 * rt_window_1x : rt_window_1x;

    int mz_start = std::max(0,
                            static_cast<int>((mz - mz_window - mz_min) / (mz_max - mz_min) * mz_bins) - 1
    );
    int mz_end = std::min(mz_bins - 1,
                          static_cast<int>((mz + mz_window - mz_min) / (mz_max - mz_min) * mz_bins) + 1
    );

    int rt_start = std::max(0,
                            static_cast<int>((rt - rt_window - rt_min) / (rt_max - rt_min) * rt_bins) - 1
    );
    int rt_end = std::min(rt_bins - 1,
                          static_cast<int>((rt + rt_window - rt_min) / (rt_max - rt_min) * rt_bins) + 1
    );

    if (rt - rt_min <= rt_window) {
      rt_start = 0;
    }
    if (rt_max - rt <= rt_window) {
      rt_end = rt_bins - 1;
    }

    for (int i = mz_start; i <= mz_end; ++i) {
      for (int j = rt_start; j <= rt_end; ++j) {
        for (int idx : cells[i][j]) {
          double mz_diff = std::abs(mz_values[idx] - mz);
          double rt_diff = std::abs(rt_values[idx] - rt);
          if (mz_diff <= mz_window && rt_diff <= rt_window) {
            result.push_back(idx);
          }
        }
      }
    }
    return result;
  }
};
double compute_noise_2d(const std::vector<int>& nearby_points,
                        const NumericVector& intensity,
                        double mz,
                        double rt,
                        double mz_window_2x,
                        double mz_window_1x,
                        double rt_window_2x,
                        double rt_window_1x,
                        const NumericVector& mz_values,
                        const NumericVector& rt_values,
                        double rt_min,
                        double rt_max) {
  double noise_sum = 0.0;
  int count = 0;

  bool is_rt_start = (rt - rt_min) <= rt_window_2x;
  bool is_rt_end = (rt_max - rt) <= rt_window_2x;

  for (int j : nearby_points) {
    double mz_diff = std::abs(mz_values[j] - mz);
    double rt_diff = std::abs(rt_values[j] - rt);

    if (is_rt_start && rt_values[j] < rt) {
      continue;
    }
    if (is_rt_end && rt_values[j] > rt) {
      continue;
    }

    bool in_noise_region = false;

    if (mz_diff >= mz_window_1x && mz_diff <= mz_window_2x) {
      in_noise_region = true;
    }

    if (rt_diff >= rt_window_1x && rt_diff <= rt_window_2x) {
      if ((is_rt_start && rt_values[j] > rt) ||
          (is_rt_end && rt_values[j] < rt) ||
          (!is_rt_start && !is_rt_end)) {
        in_noise_region = true;
      }
    }

    if (in_noise_region) {
      noise_sum += intensity[j];
      count++;
    }
  }

  return (count > 0) ? noise_sum / count : 1.0;
}

bool check_if_peak_2d(const std::vector<int>& nearby_points,
                      const NumericVector& intensity,
                      double current_intensity,
                      const NumericVector& mz_values,
                      const NumericVector& rt_values,
                      double current_mz,
                      double current_rt,
                      double mz_window,
                      double rt_window,
                      double rt_min,
                      double rt_max) {
  bool is_rt_start = (current_rt - rt_min) <= rt_window;
  bool is_rt_end = (rt_max - current_rt) <= rt_window;

  for (int j : nearby_points) {
    double mz_diff = std::abs(mz_values[j] - current_mz);
    double rt_diff = std::abs(rt_values[j] - current_rt);

    if (is_rt_start && rt_values[j] < current_rt) continue;
    if (is_rt_end && rt_values[j] > current_rt) continue;

    if (mz_diff <= mz_window && rt_diff <= rt_window && intensity[j] > current_intensity) {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
DataFrame find_2d_peaks_c(NumericVector mz,
                          NumericVector rt,
                          NumericVector intensity,
                          double mz_ppm = 5.0,
                          double rt_window = 5,
                          double snr_threshold = 3.0,
                          int mz_bins = 100,
                          int rt_bins = 100) {
  int n = mz.length();

  double mz_min = *std::min_element(mz.begin(), mz.end());
  double mz_max = *std::max_element(mz.begin(), mz.end());
  double rt_min = *std::min_element(rt.begin(), rt.end());
  double rt_max = *std::max_element(rt.begin(), rt.end());

  Grid2D grid(mz_bins, mz_min, mz_max, rt_bins, rt_min, rt_max, mz, rt);
  for (int i = 0; i < n; ++i) {
    grid.add_point(i, mz[i], rt[i]);
  }

  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&intensity](int i1, int i2) {
              return intensity[i1] > intensity[i2];
            });

  std::vector<double> result_mz;
  std::vector<double> result_rt;
  std::vector<double> result_intensity;
  std::vector<bool> checked(n, false);

  for (int i = 0; i < n; ++i) {
    int idx = indices[i];
    if (checked[idx]) continue;

    double current_mz = mz[idx];
    double current_rt = rt[idx];
    double current_intensity = intensity[idx];

    double mz_window = current_mz * mz_ppm * 1e-6;
    double mz_window_1x = mz_window;
    double mz_window_2x = 2 * mz_window;
    double rt_window_1x = rt_window;
    double rt_window_2x = 2 * rt_window;

    std::vector<int> nearby_points_2x = grid.get_nearby_points(current_mz, current_rt,
                                                               mz_window_1x, rt_window_1x,
                                                               true);  // 使用2x窗口

    double noise_mean = compute_noise_2d(nearby_points_2x, intensity, current_mz, current_rt,
                                         mz_window_2x, mz_window_1x, rt_window_2x, rt_window_1x,
                                         mz, rt, rt_min, rt_max);

    if (current_intensity < snr_threshold * noise_mean) continue;

    std::vector<int> nearby_points_1x = grid.get_nearby_points(current_mz, current_rt,
                                                               mz_window_1x, rt_window_1x,
                                                               false);  // 使用1x窗口


    bool is_peak = check_if_peak_2d(nearby_points_1x, intensity, current_intensity,
                                    mz, rt, current_mz, current_rt,
                                    mz_window_1x, rt_window_1x,
                                    rt_min, rt_max);
    for (int j : nearby_points_2x) {
      checked[j] = true;
    }

    if (is_peak) {
      result_mz.push_back(current_mz);
      result_rt.push_back(current_rt);
      result_intensity.push_back(current_intensity);
    }
  }

  return DataFrame::create(Named("mz") = result_mz,
                           Named("rt") = result_rt,
                           Named("intensity") = result_intensity);
}
