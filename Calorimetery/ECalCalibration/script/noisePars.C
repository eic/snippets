#include <DDDigi/noise/FalphaNoise.h>

double pulse(double *x, double *par) {
  return par[0] * TMath::Landau(x[0], par[2] * par[1], par[1], kTRUE);
}

void setGraphStyle(TGraph* gr, Style_t style,  Color_t color, Size_t size);

void noisePars() {

  // Signal
  const double npe_mean = 100;
  const double npe_sigma = 5;

  // To draw pulses for comparison
  constexpr double t_min = -300;
  constexpr double t_max = 300;
  constexpr double dt = 0.5;
  constexpr int n_bins = static_cast<int>((t_max - t_min) / dt);

  // Config variables
  const double gain = 55.37;
  const double sigma_offset = 3.5;
  const double sigma_analog = 10;

  const int n_samples = 7;
  const double phase = 20;
  const double time_window = 25;

  // Histograms
  TH1D* h1_adc_real = new TH1D("h1_adc_real", "", 100, 0, 450);
  TH1D* h1_adc_model = new TH1D("h1_adc_model", "", 100, 0, 450);

  // Real noise
  TRandom3 rng(0);
  const double dark_rate = 112.97e6;
  const double window_sec = (t_max - t_min) * 1.0e-9;
  const double mean_n_dark = dark_rate * window_sec;

  // Falpha noise
  const int poles = 5;
  const double variance = 1.0;
  const double alpha = 1.8;
  const double scale = 0.35;
  const double offset = 6.1;
  std::default_random_engine gen(12345);
  dd4hep::detail::FalphaNoise falpha(poles, alpha, variance);

  const int nev = 10000;

  // Example pulse shapes at a specific event
  const int evnum = 30;
  double times[n_bins] = {0};
  double amps_real[n_bins] = {0};
  double amps_model[n_bins] = {0};
  double sig_amps[n_bins] = {0};
  double noise_amps_real[n_bins] = {0};
  double noise_amps_model[n_bins] = {0};

  for (int i = 0; i < nev; i++) {
    if ((i % 500) == 0)
      cout << i << endl;

    const double sig_npe = std::max(0.0, rng.Gaus(npe_mean, npe_sigma));
    double adc_real = 0;
    double adc_model = 0;

    const int n_dark = rng.Poisson(mean_n_dark);
    std::vector<double> dark_times;
    std::vector<double> dark_amps;

    for (int i = 0; i < n_dark; i++) {
      // Dark noises are generated n_dark times randomly.
      dark_times.push_back(rng.Uniform(t_min, t_max));
      // It has 1 p.e. pulse height.
      dark_amps.push_back(1);
    }

    for (int j = 0; j < n_bins; j++) {
      double time = t_min + j * dt;
      double sig = 0;
      sig += (sig_npe * gain) * TMath::Landau(time, sigma_analog * sigma_offset,
                                              sigma_analog, kTRUE);
      double noise_real = 0;
      for (int k = 0; k < n_dark; k++) {
        noise_real +=
            (dark_amps[k] * gain) *
            TMath::Landau(time, dark_times[k] + sigma_analog * sigma_offset,
                          sigma_analog, kTRUE);
      }

      double noise_model = std::max(0.0, scale * falpha(gen) + offset);

      // Store pulse shapes for comparison
      if (i == evnum) {
        times[j] = time;
        sig_amps[j] = sig;
        noise_amps_real[j] = noise_real;
        noise_amps_model[j] = noise_model;
        amps_real[j] = sig + noise_real;
        amps_model[j] = sig + noise_model;
      }

      // CALOROC measurement
      for (int k = 0; k < n_samples; k++) {
        if (time == phase + k * time_window) {
            adc_real += sig + noise_real;
            adc_model += sig + noise_model;
        }
      }
    }
    h1_adc_real->Fill(adc_real);
    h1_adc_model->Fill(adc_model);
  }

  // For example pulses comparison
  TH2D *h2_pulse = new TH2D("h2_pulse", "", 6, -100, t_max, 6, 0, npe_mean + 4 * npe_sigma);

  TGraph *gr_sig = new TGraph(n_bins, times, sig_amps);
  setGraphStyle(gr_sig, 20, 4, 0.3);

  TGraph *gr_noise_real = new TGraph(n_bins, times, noise_amps_real);
  setGraphStyle(gr_noise_real, 20, 2, 0.3);
  TGraph *gr_noise_model = new TGraph(n_bins, times, noise_amps_model);
  setGraphStyle(gr_noise_model, 20, 2, 0.3);

  TGraph *gr_pulse_real = new TGraph(n_bins, times, amps_real);
  setGraphStyle(gr_pulse_real, 20, 1, 0.3);
  TGraph *gr_pulse_model = new TGraph(n_bins, times, amps_model);
  setGraphStyle(gr_pulse_model, 20, 1, 0.3);

  // Example pulse shapes
  TCanvas *c_pulse = new TCanvas("c_pulse", "", 1300, 1000);
  c_pulse->Divide(1, 2);

  c_pulse->cd(1);
  h2_pulse->Draw();
  gr_sig->Draw("PL");
  gr_noise_real->Draw("PL");
  gr_pulse_real->Draw("PL");

  c_pulse->cd(2);
  h2_pulse->Draw();
  gr_sig->Draw("PL");
  gr_noise_model->Draw("PL");
  gr_pulse_model->Draw("PL");

  // ADC distribution comparison
  TCanvas *c_adc = new TCanvas("c_adc", "", 1200, 600);
  c_adc->Divide(2, 1);

  c_adc->cd(1);
  h1_adc_real->Draw();

  c_adc->cd(2);
  h1_adc_model->Draw();
}

void setGraphStyle(TGraph* gr, Style_t style,  Color_t color, Size_t size){

  gr->SetMarkerStyle(style);
  gr->SetMarkerSize(size);
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
}
