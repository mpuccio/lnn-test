float weightedAverage(const std::vector<float>& w, const std::vector<float>& y) {
  float sum = 0.f;
  float total = 0.f;
  for (size_t i = 0; i < w.size() && i < y.size(); ++i) {
    sum += w[i] * y[i];
    total += w[i];
  }
  return total != 0.f ? sum / total : 0.f;
}

void valueHypertriton() {
  std::vector<float> cw{1.,2.,2.};

  std::vector<float> w0{0.5, 0.5, 1.0, 2.0};
  std::vector<float> y0{1.0128e-05, 8.5725e-06, 5.8906e-06, 1.0287e-06};
  std::vector<float> ey0{1.3441e-06, 1.1071e-06, 6.2927e-07, 1.7569e-07};
  std::vector<float> sy0{1.3294e-06, 1.5472e-06, 9.9862e-07, 3.344e-07};
  std::vector<float> csy0{8.2907e-07, 5.9867e-07, 8.092e-07, 2.9097e-07};

  std::vector<float> w1{1.0, 1.0, 2.0};
  std::vector<float> y1{2.817e-06, 2.4536e-06, 6.0039e-07};
  std::vector<float> ey1{9.4015e-07, 6.8441e-07, 2.2631e-07};
  std::vector<float> sy1{8.3386e-07, 6.5817e-07, 1.5539e-07};
  std::vector<float> csy1{6.8661e-07, 5.7153e-07, 1.033e-07};

  std::vector<float> w2{1.0, 2.0};
  std::vector<float> y2{1.9764e-06, 3.2424e-07};
  std::vector<float> ey2{3.5414e-07, 9.0408e-08};
  std::vector<float> sy2{3.8831e-07, 1.1553e-07};
  std::vector<float> csy2{2.9218e-07, 2.1032e-08};


  std::vector<float> sums{weightedAverage(w0, y0), weightedAverage(w1, y1), weightedAverage(w2, y2)};

  std::cout << "sums: " << weightedAverage(cw, sums) << std::endl;
}
