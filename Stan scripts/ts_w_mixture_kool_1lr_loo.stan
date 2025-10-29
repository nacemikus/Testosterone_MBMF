data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int<lower=1, upper=4> level1_choice[N,T];  // 1: 1 - 4
  int<lower=1, upper=2> state1[N,T];  // 1,2 
  int<lower=1, upper=2> state2[N,T];  // 1,2 
  int<lower=-4, upper=5> reward[N,T];
  
   int<lower=0, upper=1> testosterone[N];

}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;
   matrix[4, N] z;
  cholesky_factor_corr[4] L_Omega;
  
  // fixed parameters
  real beta_testo_noise;
  real beta_testo_w;
  real beta_testo_a1; 
   // real beta_testo_a2;
// 
//   real beta_testo_noise;
 
  real beta_testo_et;

  // real beta_testo;
  // real beta_ami;
  // real beta_testo;
  // real beta_ami;
  
  // Subject-level raw parameters (for Matt trick)

  
}
transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0,upper=1>[N] a1;
  vector<lower=0>[N]         beta1;
    vector<lower=0,upper=1>[N]        w;
  // vector<lower=0,upper=1>[N] a2;
  // vector<lower=0,upper=5>[N] pi;
  // vector<lower=0,upper=1>[N] w;
  vector<lower=0,upper=1>[N] lambda;
   
  matrix[4, N] r1;
  
  r1 = (diag_pre_multiply(sigma,L_Omega) * z);
  
  for (i in 1:N) {
    
      
      a1[i]     = Phi_approx(   mu_p[1] +  r1[1,i] + beta_testo_a1*testosterone[i]);
      beta1[i]  = exp( mu_p[2] +  r1[2,i] + beta_testo_noise*testosterone[i]);
   
      // a2[i]     = Phi_approx(  mu_p[5] +  r1[5,i] + beta_testo_a2*testosterone[i] );
      w[i]      = Phi_approx(mu_p[3] +  r1[3,i] + beta_testo_w*testosterone[i] );
      lambda[i] =   Phi_approx(  mu_p[4] +  r1[4,i] +  beta_testo_et*testosterone[i]);
     
    
     
  
     
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1.5);
  sigma ~ normal(0, 1); //cauchy(0,3); // 
  // fixed parameters
  beta_testo_noise ~ normal(0,1.5);

  beta_testo_a1~ normal(0,1.5);
  // beta_testo_a2~ normal(0,1.5);
  beta_testo_w ~ normal(0,1.5);
  beta_testo_et ~normal(0,1.5);

  // individual parameters
  
  to_vector(z) ~  normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);

  

  for (i in 1:N) {
    // Define values
    vector[2] v_mb;    // model-based stimulus values for level 1 (2 stimuli)
    matrix[2,2] v_mf;    // model-free stimulus values for level 1&2 (1,2--> level 1, 3-6--> level 2)
    vector[2] v_2;
    vector[2] v_hybrid;  // hybrid stimulus values for level 1 (2 stimuli)
    real level1_prob_choice2; // Initialize prob. of choosing stim 2 (0 or 1) in level 1
    matrix[2, 2] trans_prob_state1_1;
    matrix[2, 2] trans_prob_state1_2;
    int trans_prob_state1_1_changed;
    int trans_prob_state1_2_changed;
    int choice;
    int choice01;
    int choice_prev;
    vector[2] last_choice;
    
    
    // Initialize values
    trans_prob_state1_1 = rep_matrix(0.5, 2,2);
    trans_prob_state1_2 = rep_matrix(0.5, 2,2);
    v_mb  = rep_vector(0.0, 2);
    v_mf  = rep_matrix(0.0, 2,2);
    v_2 = rep_vector(0.0,2);
    v_hybrid = rep_vector(0.0, 2);
    trans_prob_state1_1_changed = 0;
    trans_prob_state1_2_changed = 0;

    for (t in 1:Tsubj[i])  {
      
      // mark that the agent has learned the state transfer
      if (state1[i,t] == 1 && trans_prob_state1_1_changed== 0) {
          trans_prob_state1_1_changed = 1;
      } else if (state1[i,t] == 2 && trans_prob_state1_2_changed== 0) {
          trans_prob_state1_2_changed = 1;
      }
      
      // compute v_mb 
      //    Qmb = Tm{s1}'*Q2;     
      if (state1[i,t]==1) {
          v_mb = trans_prob_state1_1 * v_2;
      } else if (state1[i,t]==2) {
          v_mb = trans_prob_state1_2 * v_2;
      }
      
    
      
      // compute v_hybrid
      //    Q = w*Qmb + (1-w)*Qmf(s1,:)' + st.*M(s1,:)' + respst.*R;        % mix TD and model-based value
      v_hybrid[1] = w[i] * v_mb[1] + (1-w[i]) * v_mf[state1[i,t],1];   // Q of choosing 1,3 (state 2 = 1) 
      v_hybrid[2] = w[i] * v_mb[2] + (1-w[i]) * v_mf[state1[i,t],2];   // Q of choosing 2,4 (state 2 = 2)
    
      // set the choice from 1-4 to 1-2
      choice = level1_choice[i,t];
      if(choice > 2){
        choice = choice - 2;
      }
      
      // def choice for bernoulli
      choice01 = choice -1;  
 
 //    
 
 
 
//  agent realizes transition structure:
    if (trans_prob_state1_1_changed == 1 && state1[i,t] == 1) {
         trans_prob_state1_1 = [[1, 0],[0, 1]];
    }         //  agent realizes transition structure
       
    if (trans_prob_state1_2_changed == 1 && state1[i,t] == 2) {
         trans_prob_state1_2 = [[1, 0],[0, 1]];
    } 
    
   // make the last choice sticky
  // last_choice = [0,0]';    
  //     if(!(t == 1)) {
  //       choice_prev = level1_choice[i,t-1];
  //       if(choice_prev > 2) {
  //         choice_prev = choice_prev - 2;
  //       }
  //       last_choice[choice_prev] = 1;
  //       v_hybrid[1] =  v_hybrid[1] + pi[i]*last_choice[1];
  //       v_hybrid[2] =  v_hybrid[2] + pi[i]*last_choice[2];
  //     }

     
  
       level1_prob_choice2 = inv_logit( beta1[i]*(v_hybrid[2]-v_hybrid[1]));
      choice01 ~ bernoulli(level1_prob_choice2 );  // level 1, prob. of choosing 2 in level 1
      
      // alternative model formulation
      // choice ~ categorical_logit( beta1[i]*(v_hybrid[2]-v_hybrid[1]))

      // Observe Level2 and update Level1 of the chosen option
      //    dtQ(2) = subdata.points(t) - Q2(s2);                            % prediction error (2nd choice)
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + a1[i]*(v_2[state2[i,t]] - v_mf[state1[i,t], choice] );
        
      // After observing the reward at Level 2...
      // Update Level 2 v_mf of the chosen option. Level 2--> choose one of level 2 options and observe reward
      //    Q2(s2) = Q2(s2) + lr*dtQ(2);            
      v_2[state2[i,t]] =v_2[state2[i,t]] + a1[i]*(reward[i,t] - v_2[state2[i,t]] );

      // Update Level 1 v_mf with eligibility trace
      //    Qmf(s1,a) = Qmf(s1,a) + lambda*lr*dtQ(2);   
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + lambda[i] * a1[i] * (reward[i,t] - v_2[state2[i,t]] );
      
   
    } // end of t loop
    
    
  } // end of i loop
}

generated quantities {
  // For group level parameters
  
  // real mu_sess_a1;
  // real mu_sess_beta1;
  // real mu_sess_a1;
  // real mu_sess_w;
  
  // For log likelihood calculation
  real log_lik[N*T];
  int count_trials;
  // For posterior predictive check
  real y_pred[N,T];
 
  // real y_pred_step2[N,T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N*T) {
    
    log_lik[i] = 0;
  
  }
  
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i,t] = -1;
      // y_pred_step2[i,t] = -1;
    }
  }
 
  // Generate group level parameter values
  
  { // local section, this saves time and space
  count_trials = 1;
  for (i in 1:N) {
    // Define values
       // Define values
    vector[2] v_mb;    // model-based stimulus values for level 1 (2 stimuli)
    matrix[2,2] v_mf;    // model-free stimulus values for level 1&2 (1,2--> level 1, 3-6--> level 2)
    vector[2] v_2;
    vector[2] v_hybrid;  // hybrid stimulus values for level 1 (2 stimuli)
    real level1_prob_choice2; // Initialize prob. of choosing stim 2 (0 or 1) in level 1
   
    matrix[2,2] trans_prob_state1_1;
    matrix[2,2] trans_prob_state1_2;
    int trans_prob_state1_1_changed;
    int trans_prob_state1_2_changed;
    int choice;
    int choice01;
    int choice_prev;
    vector[2] last_choice;
    // int level2_choice_01;

    // Initialize values
    trans_prob_state1_1 = rep_matrix(0.5, 2,2);
    trans_prob_state1_2 = rep_matrix(0.5, 2,2);
    v_mb  = rep_vector(0.0, 2);
    v_mf  = rep_matrix(0.0, 2,2);
    v_2 = rep_vector(0.0,2);
    v_hybrid = rep_vector(0.0, 2);
    trans_prob_state1_1_changed = 0;
    trans_prob_state1_2_changed = 0;
    
    // log_lik[i] = 0;

    for (t in 1:Tsubj[i])  {
     
        // mark that the agent has learned the state transfer
      if (state1[i,t] == 1 && trans_prob_state1_1_changed== 0) {
          trans_prob_state1_1_changed = 1;
      } else if (state1[i,t] == 2 && trans_prob_state1_2_changed== 0) {
          trans_prob_state1_2_changed = 1;
      }
      
      // compute v_mb 
      //    Qmb = Tm{s1}'*Q2;     
      if (state1[i,t]==1) {
          v_mb = trans_prob_state1_1 * v_2;
      } else if (state1[i,t]==2) {
          v_mb = trans_prob_state1_2 * v_2;
      }
      
    
      
      // compute v_hybrid
      //    Q = w*Qmb + (1-w)*Qmf(s1,:)' + st.*M(s1,:)' + respst.*R;        % mix TD and model-based value
    v_hybrid[1] = w[i] * v_mb[1] + (1-w[i]) * v_mf[state1[i,t],1];   // Q of choosing 1,3 (state 2 = 1) 
      v_hybrid[2] = w[i] * v_mb[2] + (1-w[i]) * v_mf[state1[i,t],2];   // Q of choosing 2,4 (state 2 = 2)
    
      // set the choice from 1-4 to 1-2
      choice = level1_choice[i,t];
      if(choice > 2){
        choice = choice - 2;
      }
      
        // def choice for bernoulli
      choice01 = choice -1;  
 
 //    
 
 
 
//  agent realizes transition structure:
    if (trans_prob_state1_1_changed == 1 && state1[i,t] == 1) {
         trans_prob_state1_1 = [[1, 0],[0, 1]];
    }         //  agent realizes transition structure
       
    if (trans_prob_state1_2_changed == 1 && state1[i,t] == 2) {
         trans_prob_state1_2 = [[1, 0],[0 ,1]];
    } 
    
   // make the last choice sticky
  // last_choice = [0,0]';    
  //     if(!(t == 1)) {
  //       choice_prev = level1_choice[i,t-1];
  //       if(choice_prev > 2){
  //         choice_prev = choice_prev - 2;
  //       }
  //       last_choice[choice_prev] = 1;
  //       v_hybrid[1] =  v_hybrid[1] + pi[i]*last_choice[1];
  //       v_hybrid[2] =  v_hybrid[2] + pi[i]*last_choice[2];
  //     }
      
    
         level1_prob_choice2 = inv_logit( beta1[i]*(v_hybrid[2]-v_hybrid[1]));
  
      log_lik[count_trials] = log_lik[count_trials] + bernoulli_lpmf( choice01 | level1_prob_choice2 );
      
      if (is_nan(log_lik[count_trials])) log_lik[count_trials] =0;
      
      count_trials = count_trials +1;
      
      
       // generate posterior prediction for current trial
      y_pred[i,t] = bernoulli_rng(level1_prob_choice2);

      // alternative model formulation
      // choice ~ categorical_logit( beta1[i]*(v_hybrid[2]-v_hybrid[1]))

      // Observe Level2 and update Level1 of the chosen option
      //    dtQ(2) = subdata.points(t) - Q2(s2);                            % prediction error (2nd choice)
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + a1[i]*(v_2[state2[i,t]] - v_mf[state1[i,t], choice] );
        
      // After observing the reward at Level 2...
      // Update Level 2 v_mf of the chosen option. Level 2--> choose one of level 2 options and observe reward
      //    Q2(s2) = Q2(s2) + lr*dtQ(2);            
      v_2[state2[i,t]] =v_2[state2[i,t]] + a1[i]*(reward[i,t] - v_2[state2[i,t]] );

      // Update Level 1 v_mf with eligibility trace
       //    Qmf(s1,a) = Qmf(s1,a) + lambda*lr*dtQ(2);   
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + lambda[i] * a1[i] * (reward[i,t] - v_2[state2[i,t]] );
        } // end of t loop
        
       
    } // end of i loop
   } // end local
}

