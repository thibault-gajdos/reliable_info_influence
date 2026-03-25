functions {
  vector clamp_vector(vector x, real lo, real hi) {
    vector[num_elements(x)] out;
    for (i in 1:num_elements(x)) out[i] = fmin(fmax(x[i], lo), hi);
    return out;
  }

  real clamp_real(real x, real lo, real hi) {
    return fmin(fmax(x, lo), hi);
  }

  // ===== Simple mapping: m(p) = alpha * logit(p) + beta =====
  real map_linear_logodds(real p, real alpha, real beta) {
    return alpha * logit(p) + beta;
  }

  // ===== Sequential weight (w6 = 1 fixed) =====
  real get_weight(int s, vector w5) {
    if (s <= 5) return w5[s];
    return 1.0;
  }

  // Worker function for reduce_sum
  real partial_sum(array[] int slice_indices,
                   int start, int end,
                   vector mu_pr, vector sigma_pr,
                   array[] int Tsubj,
                   array[,] int sample,
                   array[,,] int color,
                   array[,,] real proba,
                   array[,] int choice,
                   array[,,] int influence_sample,
                   array[,] real logit_influence,
                   matrix param_raw) {
    real lp = 0;

    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];

      // Subject-level parameters:
      // 1-5: w1..w5 in [0, 2]
      // 6:   alpha in [0, 10]
      // 7:   beta (unbounded)
      // 8:   a_infl (unbounded)
      // 9:   b_infl (unbounded)
      // 10:  sigma_infl > 0

      vector[5] w5;
      for (j in 1:5)
        w5[j] = 2 * Phi_approx(mu_pr[j] + sigma_pr[j] * param_raw[n, j]);

      real alpha      = 10 * Phi_approx(mu_pr[6] + sigma_pr[6] * param_raw[n, 6]);
      real beta       = mu_pr[7] + sigma_pr[7] * param_raw[n, 7];
      real a_infl     = mu_pr[8] + sigma_pr[8] * param_raw[n, 8];
      real b_infl     = mu_pr[9] + sigma_pr[9] * param_raw[n, 9];
      real sigma_infl = exp(clamp_real(mu_pr[10] + sigma_pr[10] * param_raw[n, 10], -20, 20));

      for (t in 1:Tsubj[n]) {
        int S = sample[n, t];
        if (S < 1) continue;

        int ch = choice[n, t];
        if (ch < 1 || ch > 2) continue;

        vector[2] evidence = rep_vector(0.0, 2);
        real sub_evidence = 0.0;

        for (s in 1:S) {
          real p = proba[n, t, s];
          int c  = color[n, t, s];

          if (c < 1 || c > 2) continue;
          if (p <= 0 || p >= 1) continue;
          if (is_nan(p) || is_inf(p)) continue;

          real m  = map_linear_logodds(p, alpha, beta);
          real ws = get_weight(s, w5);

          evidence[c] += ws * m;

          // Sub-evidence for influence-matching samples
          if (influence_sample[n, t, s] == 1) {
            real d_color = (c == 1) ? 1.0 : -1.0;
            sub_evidence += ws * m * d_color;
          }
        }

        // --- Choice likelihood ---
        vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
        lp += bernoulli_logit_lpmf(ch == 1 | evidence_safe[1] - evidence_safe[2]);

        // --- Influence report likelihood ---
        // Flip both sub and total evidence relative to choice
        real ev_diff = evidence_safe[1] - evidence_safe[2];
        real sub_ev_choice   = (ch == 1) ? sub_evidence : -sub_evidence;
        real total_ev_choice = (ch == 1) ? ev_diff : -ev_diff;

        // Relative influence: ratio of sub to total, guarded
        real denom = fmax(total_ev_choice, 1e-4);
        real rel_influence = clamp_real(sub_ev_choice / denom, -100, 100);

        real z_obs = logit_influence[n, t];
        if (z_obs > -90) {
          real mu_infl = clamp_real(a_infl * rel_influence + b_infl, -100, 100);
          lp += normal_lpdf(z_obs | mu_infl, sigma_infl);
        }
      }
    }

    return lp;
  }

  // ===== Helper: compute evidence + sub_evidence for a trial =====
  vector compute_evidence_and_sub(int sample_size,
                                  array[] int color_data,
                                  array[] real proba_data,
                                  array[] int influence_sample_data,
                                  vector w5,
                                  real alpha, real beta) {
    // out[1] = evidence[1], out[2] = evidence[2], out[3] = sub_evidence
    vector[3] out = rep_vector(0.0, 3);

    for (s in 1:sample_size) {
      real p = proba_data[s];
      int c  = color_data[s];

      if (c < 1 || c > 2) continue;
      if (p <= 0 || p >= 1) continue;
      if (is_nan(p) || is_inf(p)) continue;

      real m  = map_linear_logodds(p, alpha, beta);
      real ws = get_weight(s, w5);

      out[c] += ws * m;

      if (influence_sample_data[s] == 1) {
        real d_color = (c == 1) ? 1.0 : -1.0;
        out[3] += ws * m * d_color;
      }
    }
    return out;
  }
}

data {
  int<lower=1> N;
  int<lower=1> T_max;
  int<lower=1> I_max;

  array[N] int<lower=1> Tsubj;
  array[N, T_max] int<lower=-1> sample;
  array[N, T_max, I_max] int<lower=-1, upper=2> color;
  array[N, T_max, I_max] real<lower=-1, upper=1> proba;
  array[N, T_max] int<lower=-1, upper=2> choice;

  // Influence data
  array[N, T_max, I_max] int<lower=-1, upper=1> influence_sample;
  array[N, T_max] real logit_influence;  // pre-transformed, use -99 for missing

  int<lower=5> grainsize;
}

parameters {
  vector[10] mu_pr;
  vector<lower=0>[10] sigma_pr;
  matrix[N, 10] param_raw;
}

model {
  mu_pr ~ std_normal();
  sigma_pr ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, choice,
                       influence_sample, logit_influence,
                       param_raw);
}

generated quantities {
  matrix[N, 10] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max);
  array[N, T_max] real pred_proba = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);
  array[N, T_max] real rel_infl_out = rep_array(-99.0, N, T_max);
  array[N, T_max] real sub_ev_out = rep_array(-99.0, N, T_max);
  array[N, T_max] real pred_logit_influence = rep_array(-99.0, N, T_max);
  array[N, T_max] real pred_influence = rep_array(-99.0, N, T_max);

  vector[sum(Tsubj)] log_lik;
  vector[sum(Tsubj)] log_lik_choice;
  vector[sum(Tsubj)] log_lik_influence;

  // Interpret group-level means
  real mu_w1    = 2 * Phi_approx(mu_pr[1]);
  real mu_w2    = 2 * Phi_approx(mu_pr[2]);
  real mu_w3    = 2 * Phi_approx(mu_pr[3]);
  real mu_w4    = 2 * Phi_approx(mu_pr[4]);
  real mu_w5    = 2 * Phi_approx(mu_pr[5]);
  real mu_alpha = 10 * Phi_approx(mu_pr[6]);
  real mu_beta  = mu_pr[7];
  real mu_a     = mu_pr[8];
  real mu_b     = mu_pr[9];
  real mu_sigma_infl = exp(mu_pr[10]);

  int k = 0;
  for (n in 1:N) {
    // Transform subject-level parameters
    vector[5] w5;
    for (j in 1:5) {
      w5[j] = 2 * Phi_approx(mu_pr[j] + sigma_pr[j] * param_raw[n, j]);
      params[n, j] = w5[j];
    }
    params[n, 6]  = 10 * Phi_approx(mu_pr[6] + sigma_pr[6] * param_raw[n, 6]);
    params[n, 7]  = mu_pr[7] + sigma_pr[7] * param_raw[n, 7];
    params[n, 8]  = mu_pr[8] + sigma_pr[8] * param_raw[n, 8];
    params[n, 9]  = mu_pr[9] + sigma_pr[9] * param_raw[n, 9];
    params[n, 10] = exp(clamp_real(mu_pr[10] + sigma_pr[10] * param_raw[n, 10], -20, 20));

    for (t in 1:Tsubj[n]) {
      k += 1;

      if (sample[n, t] < 1) {
        y_pred[n, t] = -1.0;
        pred_proba[n, t] = -1.0;
        diff_evidence[n, t] = -1.0;
        log_lik[k] = 0;
        log_lik_choice[k] = 0;
        log_lik_influence[k] = 0;
        continue;
      }

      // Pack trial arrays
      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      array[I_max] int infl_sample_trial;
      for (ii in 1:I_max) {
        color_trial[ii] = color[n, t, ii];
        proba_trial[ii] = proba[n, t, ii];
        infl_sample_trial[ii] = influence_sample[n, t, ii];
      }

      vector[3] ev = compute_evidence_and_sub(sample[n, t],
                                               color_trial, proba_trial,
                                               infl_sample_trial,
                                               w5,
                                               params[n, 6], params[n, 7]);

      // --- Choice ---
      vector[2] evidence;
      evidence[1] = ev[1];
      evidence[2] = ev[2];
      vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
      real logit_diff = evidence_safe[1] - evidence_safe[2];
      real p_blue = inv_logit(logit_diff);

      y_pred[n, t] = bernoulli_rng(p_blue) + 1;  // +1 to keep {1,2} encoding
      pred_proba[n, t] = p_blue;
      diff_evidence[n, t] = logit_diff;

      real ll_choice = bernoulli_logit_lpmf(choice[n, t] == 1 | logit_diff);

      // --- Influence report ---
      int ch = choice[n, t];
      real sub_ev_choice   = (ch == 1) ? ev[3] : -ev[3];
      real total_ev_choice = (ch == 1) ? logit_diff : -logit_diff;

      real denom = fmax(total_ev_choice, 1e-4);
      real rel_influence = clamp_real(sub_ev_choice / denom, -100, 100);

      sub_ev_out[n, t] = sub_ev_choice;
      rel_infl_out[n, t] = rel_influence;

      real ll_infl = 0;
      real z_obs = logit_influence[n, t];

      if (z_obs > -90) {
        real mu_infl = clamp_real(params[n, 8] * rel_influence + params[n, 9], -100, 100);

        ll_infl = normal_lpdf(z_obs | mu_infl, params[n, 10]);

        // Posterior predictive
        real z_pred = normal_rng(mu_infl, params[n, 10]);
        pred_logit_influence[n, t] = z_pred;
        pred_influence[n, t] = inv_logit(z_pred) - 0.5;
      }

      log_lik_choice[k] = ll_choice;
      log_lik_influence[k] = ll_infl;
      log_lik[k] = ll_choice + ll_infl;
    }
  }
}
