functions {
  vector clamp_vector(vector x, real lo, real hi) {
    vector[num_elements(x)] out;
    for (i in 1:num_elements(x)) out[i] = fmin(fmax(x[i], lo), hi);
    return out;
  }

  // ===== Simple mapping: m(p) = alpha * logit(p) + beta =====
  real map_linear_logodds(real p, real alpha, real beta) {
    return alpha * logit(p) + beta;
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
                   matrix param_raw) {
    real lp = 0;

    for (i in 1:size(slice_indices)) {
      int n = slice_indices[i];

      // Subject-level parameters (still 5 to keep structure unchanged):
      // 1 = lambda_recency in (0,1)
      // 2 = alpha in [0,2]
      // 3 = beta gaussian
      vector[4] params;

      params[1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw[n, 1]);        // lambda_recency
      params[2] = 10 * Phi_approx(mu_pr[2] + sigma_pr[2] * param_raw[n, 2]);    // alpha in [0,2]
      params[3] = mu_pr[3] + sigma_pr[3] * param_raw[n, 3]; 

      for (t in 1:Tsubj[n]) {
        int S = sample[n, t];
        if (S < 1) continue;

        vector[2] evidence = rep_vector(0.0, 2);

        for (s in 1:S) {
          real p = proba[n, t, s];
          int c  = color[n, t, s];

          if (c < 1 || c > 2) continue;
          if (p <= 0 || p >= 1) continue;
          if (is_nan(p) || is_inf(p)) continue;

          // linear log-odds contribution
          real m = map_linear_logodds(p, params[2], params[3]);

          evidence[c] += exp(params[1] * (s - S)) * m;
        }

        vector[2] evidence_safe = clamp_vector(evidence, -100, 100);

        // theta fixed to 1: softmax(evidence_safe)
        lp += bernoulli_logit_lpmf(choice[n, t] == 1 | evidence_safe[1] - evidence_safe[2]);
      }
    }

    return lp;
  }

  vector compute_evidence(int sample_size,
                          array[] int color_data,
                          array[] real proba_data,
                          real lambda_recency,
                          real alpha, real beta
) {    
    vector[2] evidence = rep_vector(0.0, 2);

    for (s in 1:sample_size) {
      real p = proba_data[s];
      int c  = color_data[s];

      if (c < 1 || c > 2) continue;
      if (p <= 0 || p >= 1) continue;
      if (is_nan(p) || is_inf(p)) continue;

      real m = map_linear_logodds(p, alpha, beta);
      evidence[c] += exp(lambda_recency * (s - sample_size)) * m;
    }
    return evidence;
  }

  real compute_log_lik(int sample_size,
                       array[] int color_data,
                       array[] real proba_data,
		       int choice_obs,
                       real lambda_recency,
                       real alpha, real beta) {    
    vector[2] evidence = compute_evidence(sample_size, color_data, proba_data,
                                         lambda_recency, alpha, beta);
    vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
   return bernoulli_logit_lpmf(choice_obs == 1 | evidence_safe[1] - evidence_safe[2]);
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

  int<lower=5> grainsize;
}

parameters {
  vector[3] mu_pr;
  vector<lower=0>[3] sigma_pr;
  matrix[N, 3] param_raw;
}

model {
  mu_pr ~ std_normal();
  sigma_pr ~ normal(0, 0.5);
  to_vector(param_raw) ~ std_normal();

  array[N] int indices;
  for (n in 1:N) indices[n] = n;

  target += reduce_sum(partial_sum, indices, grainsize,
                       mu_pr, sigma_pr,
                       Tsubj, sample, color, proba, choice, param_raw);
}

generated quantities {
  matrix[N, 3] params;
  array[N, T_max] real y_pred = rep_array(-1.0, N, T_max);
  array[N, T_max] real pred_proba = rep_array(-1.0, N, T_max);
  array[N, T_max] real diff_evidence = rep_array(-1.0, N, T_max);
  vector[sum(Tsubj)] log_lik;

  // interpret mu's (consistent with transforms)
  real mu_lambda = Phi_approx(mu_pr[1]);
  real mu_alpha  = 10 * Phi_approx(mu_pr[2]);
  real mu_beta   = mu_pr[3];

  int k = 0;
  for (n in 1:N) {
    params[n, 1] = Phi_approx(mu_pr[1] + sigma_pr[1] * param_raw[n, 1]);          // lambda_recency
    params[n, 2] = 10 * Phi_approx(mu_pr[2] + sigma_pr[2] * param_raw[n, 2]);      // alpha
    params[n, 3] = mu_pr[3] + sigma_pr[3] * param_raw[n, 3]; // beta

    for (t in 1:Tsubj[n]) {
      k += 1;

      if (sample[n, t] < 1) {
        y_pred[n, t] = -1.0;
        pred_proba[n, t] = -1.0;
        diff_evidence[n, t] = -1.0;
        log_lik[k] = 0;
        continue;
      }

      array[I_max] int color_trial;
      array[I_max] real proba_trial;
      for (i in 1:I_max) {
        color_trial[i] = color[n, t, i];
        proba_trial[i] = proba[n, t, i];
      }

      vector[2] evidence = compute_evidence(sample[n, t], color_trial, proba_trial,
                                            params[n, 1],
                                            params[n, 2], params[n, 3]);

      vector[2] evidence_safe = clamp_vector(evidence, -100, 100);
      real logit_diff = evidence_safe[1] - evidence_safe[2];
      real p_blue = inv_logit(logit_diff);
      y_pred[n, t] = bernoulli_rng(p_blue) + 1;  // +1 to keep {1,2} encoding
      pred_proba[n, t] = p_blue;
      diff_evidence[n, t] = evidence_safe[1] - evidence_safe[2];

      log_lik[k] = bernoulli_logit_lpmf(choice[n, t] == 1 | logit_diff);
    }
  }
}
