#define TMB_LIB_INIT R_init_RBassess
#include <TMB.hpp>

// Calculate abundance-at-age and length
template<class Type>
Type get_N(int a, Type Lmid, Type R, Type F, Type M, Type Linf, Type gamma_angler, Type SL50_angler,
           vector<Type> Len_age) {
  Type ans = R;
  if(a>0) {
    Type F_accumulator = 0;
    Type M_accumulator = 0;
    Type rel_length = Lmid/Len_age(a);
    for(int a_prime=0;a_prime<a;a_prime++) {
      F_accumulator += pow(1 + exp(-gamma_angler * (rel_length * Len_age(a_prime) - SL50_angler)), -1);
      M_accumulator += M;
    }
    F_accumulator *= F;
    ans *= exp(-M_accumulator - F_accumulator);
  }
  return CppAD::CondExpGt(ans, Type(1e-8), ans, Type(1e-8));
}

// Baranov equation
template<class Type>
Type get_C(int len, Type N, Type F, Type M, vector<Type> sel) {
  Type Z = sel(len) * F + M;
  Type answer = sel(len) * F * N * (1 - exp(-Z));
  return answer/Z;
}

// Calculate selectivity at length
template<class Type>
vector<Type> calculate_logistic_selectivity(vector<Type> Lmid, Type gamma, Type SL50) {
  int n_bins = Lmid.size();
  vector<Type> sel(n_bins);
  for(int len=0;len<n_bins;len++) sel(len) = pow(1 + exp(-gamma * (Lmid(len) - SL50)), -1);
  Type max_sel = max(sel);
  for(int len=0;len<n_bins;len++) sel(len) /= max_sel;
  return sel;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Len_data);
  DATA_MATRIX(Len_age_data);
  DATA_VECTOR(age_adjust);
  DATA_VECTOR(Lmid);
  DATA_SCALAR(Lbin_width);

  DATA_VECTOR(stocking_density);
  DATA_VECTOR(L_stock);
  DATA_SCALAR(bag_limit);
  DATA_SCALAR(release_mortality);
  DATA_SCALAR(p_vrel);
  DATA_SCALAR(init_p_harvest);
  DATA_INTEGER(nit_harvest);

  DATA_INTEGER(use_likelihood);
  DATA_INTEGER(use_priors);

  DATA_VECTOR(prior_Linf);
  DATA_VECTOR(prior_K);
  DATA_VECTOR(prior_CV_Len);
  DATA_VECTOR(prior_M);

  DATA_VECTOR(prior_Effort);
  DATA_VECTOR(prior_q);

  DATA_VECTOR(prior_GN_SL50);
  DATA_VECTOR(prior_GN_gamma);
  DATA_VECTOR(prior_angler_SL50);
  DATA_VECTOR(prior_angler_gamma);

  PARAMETER(Linf);
  PARAMETER(K);
  PARAMETER(CV_Len);
  PARAMETER(M);

  PARAMETER(Effort);
  PARAMETER(q);

  PARAMETER(GN_SL50);
  PARAMETER(GN_gamma);
  PARAMETER(angler_SL50);
  PARAMETER(angler_gamma);

  // Calculate length-at-age
  int n_age = age_adjust.size();
  vector<Type> Len_age(n_age);
  vector<Type> sd_Len_age(n_age);

  for(int a=0;a<n_age;a++) {
    Len_age(a) = L_stock(a) * exp(-K * age_adjust(a)) + Linf * (1 - exp(-K * age_adjust(a)));
    sd_Len_age(a) = CV_Len * Len_age(a);
  }

  // Selectivity (gillnet and angler)
  int n_bins = Lmid.size();
  vector<Type> sel_GN(n_bins);
  vector<Type> sel_angler(n_bins);

  sel_GN = calculate_logistic_selectivity(Lmid, GN_gamma, GN_SL50);
  sel_angler = calculate_logistic_selectivity(Lmid, angler_gamma, angler_SL50);

  // Population model
  Type p_harvest;
  Type CPUE;
  matrix<Type> N(n_age,n_bins);
  matrix<Type> Cat(n_age,n_bins);
  matrix<Type> survey_N(n_age,n_bins);

  vector<Type> survey_NL(n_bins);
  survey_NL.setZero();

  // Initial value for F and iterate over j to solve for F
  vector<Type> F_vec(nit_harvest);
  F_vec(0) = q * Effort * (init_p_harvest + (1 - init_p_harvest) * release_mortality);
  for(int j=1;j<nit_harvest;j++) {

    for(int a=0;a<n_age;a++) {
      for(int len=0;len<n_bins;len++) {
        // Length-at-age probability
        Type zn = pnorm(Lmid(len) + 0.5 * Lbin_width, Len_age(a), sd_Len_age(a));
        zn -= pnorm(Lmid(len) - 0.5 * Lbin_width, Len_age(a), sd_Len_age(a));

        // Calculate abundance at age and length
        N(a,len) = get_N(a, Lmid(len), zn * stocking_density(a), F_vec(j-1), M, Linf, angler_gamma, angler_SL50, Len_age);

        // Baranov equation
        Cat(a,len) = get_C(len, N(a,len), F_vec(j-1), M, sel_angler);

        // Surveyed abundance
        survey_N(a,len) = sel_GN(len) * N(a,len);
        survey_NL(len) += survey_N(a,len);
      }
    }
    CPUE = Cat.sum()/Effort;

    // Calculate retention rate from CPUE subject to bag limit and voluntary release
    p_harvest = 0;
    for(int i=0;i<=20;i++) {
      Type x = i;
      Type max_x = CppAD::CondExpGt(x * (1 - p_vrel), bag_limit, bag_limit, x * (1 - p_vrel));
      p_harvest += max_x * exp(-CPUE) * pow(CPUE, x)/exp(lfactorial(x));
    }
    p_harvest /= CPUE;

    // New value for F in the j-th iteration
    F_vec(j) = q * Effort * (p_harvest + (1 - p_harvest) * release_mortality);
  }

  Type F = F_vec(nit_harvest-1);
  Type F_retain = q * Effort * p_harvest;
  Type F_release = F - F_retain;

  // Negative log prior density
  Type neg_log_prior = 0;
  if(use_priors) {
    neg_log_prior -= dnorm(log(Linf), prior_Linf(0), prior_Linf(1), true);
    neg_log_prior -= dnorm(log(K), prior_K(0), prior_K(1), true);
    neg_log_prior -= dnorm(log(CV_Len), prior_CV_Len(0), prior_CV_Len(1), true);
    neg_log_prior -= dnorm(log(M), prior_M(0), prior_M(1), true);

    neg_log_prior -= dnorm(log(Effort), prior_Effort(0), prior_Effort(1), true);
    neg_log_prior -= dnorm(log(q), prior_q(0), prior_q(1), true);

    neg_log_prior -= dnorm(log(GN_SL50), prior_GN_SL50(0), prior_GN_SL50(1), true);
    neg_log_prior -= dnorm(log(GN_gamma), prior_GN_gamma(0), prior_GN_gamma(1), true);
    neg_log_prior -= dnorm(log(angler_SL50), prior_angler_SL50(0), prior_angler_SL50(1), true);
    neg_log_prior -= dnorm(log(angler_gamma), prior_angler_gamma(0), prior_angler_gamma(1), true);
  }

  // Negative log likelihood
  vector<Type> neg_log_like(2);
  neg_log_like.setZero();
  if(use_likelihood) {
    Type sum_survey_N = survey_N.sum();
    Type sum_survey_NL = survey_NL.sum();
    for(int len=0;len<n_bins;len++) {
      for(int a=0;a<n_age;a++) {
        if(!R_IsNA(asDouble(Len_age_data(a,len))) && Len_age_data(a,len) > 0) {
          neg_log_like(0) -= Len_age_data(a,len) * log(survey_N(a,len)/sum_survey_N);
        }
      }
      if(!R_IsNA(asDouble(Len_data(len))) && Len_data(len) > 0) {
        neg_log_like(1) -= Len_data(len) * log(survey_NL(len)/sum_survey_NL);
      }
    }
  }

  // Negative log posterior
  Type neg_log_posterior = neg_log_prior + neg_log_like.sum();

  // Report model quantities to R
  ADREPORT(F);
  ADREPORT(p_harvest);
  ADREPORT(CPUE);

  REPORT(F_vec);
  REPORT(F_retain);
  REPORT(F_release);

  REPORT(Len_age);
  REPORT(Linf);
  REPORT(K);
  REPORT(CV_Len);
  REPORT(sd_Len_age);

  REPORT(GN_gamma);
  REPORT(GN_SL50);
  REPORT(angler_gamma);
  REPORT(angler_SL50);

  REPORT(F);
  REPORT(Effort);
  REPORT(q);

  REPORT(CPUE);
  REPORT(p_harvest);

  REPORT(N);
  REPORT(Cat);
  REPORT(survey_N);
  REPORT(survey_NL);

  REPORT(sel_GN);
  REPORT(sel_angler);

  REPORT(neg_log_prior);
  REPORT(neg_log_like);
  REPORT(neg_log_posterior);

  return neg_log_posterior;
}



