data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int<lower=1, upper=4> level1_choice[N,T];  // 1: 1 - 4
  int<lower=1, upper=2> state1[N,T];  // 1,2 
  int<lower=1, upper=2> state2[N,T];  // 1,2 
  int<lower=1, upper=4> stim_left[N,T];  // 1-4
  int<lower=-4, upper=5> reward[N,T];
  
  // int<lower=1, upper=T> Tsubj_t0[N];
  // int<lower=-1, upper=4> level1_choice_t0[N,T];  // 1: 1 - 4
  // int<lower=-1, upper=2> state1_t0[N,T];  // 1,2 
  // int<lower=-1, upper=2> state2_t0[N,T];  // 1,2 
  //  int<lower=1, upper=4> stim_left_t0[N,T];  // 1-4
  // 
  // int<lower=-4, upper=5> reward_t0[N,T];
  
 
  int<lower=0, upper=1> testosterone[N];
  real dat1[N];
  real comt[N];
  real cag[N];
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[3] mu_p;
  vector<lower=0>[3] sigma;
   matrix[3, N] z;
  cholesky_factor_corr[3] L_Omega;
  // fixed parameters
  real beta_testosterone_w;
 
  real beta_testosterone_g;

  real beta_testosterone_noise;
  
  real beta_comt_w;
  real beta_comt_noise;
  real beta_comt_g;
  
  real beta_dat1_w;
  real beta_dat1_noise;
  real beta_dat1_g;
  
  real beta_cag_w;
  real beta_cag_noise;
  real beta_cag_g;
  
  real beta_comt_testosterone_w;
  real beta_comt_testosterone_noise;
  real beta_comt_testosterone_g;
  
    
  real beta_dat1_testosterone_w;
  real beta_dat1_testosterone_noise;
  real beta_dat1_testosterone_g;
    
  real beta_cag_testosterone_w;
  real beta_cag_testosterone_noise;
  real beta_cag_testosterone_g;
  

}
transformed parameters {
  // Transform subject-level raw parameters
 
  vector<lower=0,upper=1>[N] g;
  vector<lower=0>[N] noise;

  vector<lower=0,upper=1>[N] w;
  
  matrix[3, N] r1;
  
  r1 = (diag_pre_multiply(sigma,L_Omega) * z);
  
  for (i in 1:N) {
      
      
      
    w[i]     = Phi_approx(mu_p[1]+ r1[1,i] + (beta_comt_testosterone_w*comt[i] + 
    beta_dat1_testosterone_w*dat1[i] + beta_cag_testosterone_w*cag[i] + beta_testosterone_w)*testosterone[i] +
    beta_comt_w*comt[i] + beta_dat1_w*dat1[i] + beta_cag_w*cag[i]);
    
    noise[i]     =  exp( mu_p[2]+ r1[2,i]  + (beta_comt_testosterone_noise*comt[i] + 
    beta_dat1_testosterone_noise*dat1[i] + beta_cag_testosterone_noise*cag[i] + beta_testosterone_noise)*testosterone[i] +
    beta_comt_noise*comt[i] + beta_dat1_noise*dat1[i] + beta_cag_noise*cag[i]);
    
    g[i] =Phi_approx( mu_p[3]+ r1[3,i]  + (beta_comt_testosterone_g*comt[i] + 
    beta_dat1_testosterone_g*dat1[i] + beta_cag_testosterone_g*cag[i] + beta_testosterone_g)*testosterone[i] +
    beta_comt_g*comt[i] + beta_dat1_g*dat1[i] + beta_cag_g*cag[i]);
     
    
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1.5);
 
  sigma ~ normal(0,1);
 
  // fixed parameters
  beta_testosterone_w ~ normal(0,1.5);
 
  beta_testosterone_g ~ normal(0,1.5);
 
  beta_testosterone_noise ~ normal(0,1.5);
  beta_comt_w ~ normal(0,1.5);
  beta_comt_noise ~ normal(0,1.5);
  beta_comt_g ~ normal(0,1.5);
  
  beta_dat1_w ~ normal(0,1.5);
  beta_dat1_noise ~ normal(0,1.5);
  beta_dat1_g ~ normal(0,1.5);
  
  beta_cag_w ~ normal(0,1.5);
  beta_cag_noise ~ normal(0,1.5);
  beta_cag_g ~ normal(0,1.5);
  
  beta_comt_testosterone_w ~ normal(0,1.5);
  beta_comt_testosterone_noise ~ normal(0,1.5);
  beta_comt_testosterone_g ~ normal(0,1.5);
  
    
  beta_dat1_testosterone_w ~ normal(0,1.5);
  beta_dat1_testosterone_noise ~ normal(0,1.5);
  beta_dat1_testosterone_g ~ normal(0,1.5);
    
  beta_cag_testosterone_w ~ normal(0,1.5);
  beta_cag_testosterone_noise ~ normal(0,1.5);
  beta_cag_testosterone_g ~ normal(0,1.5);
  
  
  to_vector(z) ~  normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  // individual parameters
 

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
    int action_right;
    int action_left;
    int pressed_left_prev;
    
    real a1_par;
    real a1_g_par;
    real a1_l_par;
    real a2_par;
    real g_par;
    real st_par;
    real w_par;
    real delta;
    
    w_par = w[i];
    a1_g_par = 1;
    a1_l_par = 1;
    a2_par =1;
    g_par = g[i];
    st_par = 0;
    
    
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
      v_hybrid[1] = w_par * v_mb[1] + (1-w_par) * v_mf[state1[i,t],1];   // Q of choosing 1,3 (state 2 = 1) 
      v_hybrid[2] = w_par * v_mb[2] + (1-w_par) * v_mf[state1[i,t],2];   // Q of choosing 2,4 (state 2 = 2)
    
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
    
  // make the last response sticky
  
      if(!(t == 1)) {
        pressed_left_prev = stim_left[i,t-1] == level1_choice[i,t-1];
        action_left = stim_left[i,t];
        if (action_left > 2) action_left = action_left- 2;
        
        action_right = 3 - action_left; 
        
        if(pressed_left_prev == 1) {
          v_hybrid[action_left] =  v_hybrid[action_left] + st_par;
        } else {
           v_hybrid[action_right] =  v_hybrid[action_right] + st_par;
        }
       
      }

      
      level1_prob_choice2 = inv_logit( noise[i]*(v_hybrid[2]-v_hybrid[1]));
       // level1_prob_choice2 = noise[i]*(v_hybrid[2]-v_hybrid[1]);
  
      choice01 ~ bernoulli(level1_prob_choice2 );  // level 1, prob. of choosing 2 in level 1
      
      // alternative model formulation
      if (is_nan(sum(noise[i]*v_hybrid))) {
        print("beta is ", noise[i], ", vh[1] is ", v_hybrid[1], ", vh[2] is ", v_hybrid[1] );
      }
      // choice ~ categorical_logit(noise[i]*v_hybrid);
      

      // Observe Level2 and update Level1 of the chosen option
      //    dtQ(2) = subdata.points(t) - Q2(s2);                            % prediction error (2nd choice)
      delta = reward[i,t] - v_mf[state1[i,t], choice];
      if (delta >0) {
        a1_par = a1_g_par;
      } else {
        a1_par = a1_l_par;
      }
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + a1_par*delta;
        
      // After observing the reward at Level 2...
      // Update Level 2 v_mf of the chosen option. Level 2--> choose one of level 2 options and observe reward
      //    Q2(s2) = Q2(s2) + lr*dtQ(2);            
      v_2[state2[i,t]] =v_2[state2[i,t]] + a2_par*(reward[i,t] - v_2[state2[i,t]] );

      // Update Level 1 v_mf with eligibility trace
      //    Qmf(s1,a) = Qmf(s1,a) + lambda*lr*dtQ(2);   
      // v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + lambda[i] * a1[i] * (reward[i,t] - v_2[state2[i,t]] );
      // forget others
      
      v_mf[state1[i,t], 3-choice] = (1-g_par)*v_mf[state1[i,t], 3-choice];
      v_mf[3-state1[i,t], choice] = (1-g_par)*v_mf[3-state1[i,t], choice];
      v_mf[3-state1[i,t], 3-choice] = (1-g_par )*v_mf[3-state1[i,t], 3-choice];
      
   
    } // end of t loop
    
   
    
  } // end of i loop
}

generated quantities {
 
  
  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N,T];
  // real y_pred_t0[N,T];
  // real y_pred_step2[N,T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i,t] = -1;
      // y_pred_step2[i,t] = -1;
    }
  }
  
  { // local section, this saves time and space
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
    int action_right;
    int action_left;
    int pressed_left_prev;
    // int level2_choice_01;
    
    real delta;
    real a1_par;
      real a1_g_par;
      real a1_l_par;
    real a2_par;
    real g_par;
    real st_par;
    real w_par;
    
    w_par=w[i];
    a1_g_par = 1;
    a1_l_par = 1;
    a2_par = 1;
    g_par = g[i];
    st_par = 0;
    // Initialize values
    trans_prob_state1_1 = rep_matrix(0.5, 2,2);
    trans_prob_state1_2 = rep_matrix(0.5, 2,2);
    v_mb  = rep_vector(0.0, 2);
    v_mf  = rep_matrix(0.0, 2,2);
    v_2 = rep_vector(0.0,2);
    v_hybrid = rep_vector(0.0, 2);
    trans_prob_state1_1_changed = 0;
    trans_prob_state1_2_changed = 0;
    
    log_lik[i] = 0;

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
      v_hybrid[1] = w_par * v_mb[1] + (1-w_par) * v_mf[state1[i,t],1];   // Q of choosing 1,3 (state 2 = 1) 
      v_hybrid[2] = w_par * v_mb[2] + (1-w_par) * v_mf[state1[i,t],2];   // Q of choosing 2,4 (state 2 = 2)
    
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
    if(!(t == 1)) {
        pressed_left_prev = stim_left[i,t-1] == level1_choice[i,t-1];
        action_left = stim_left[i,t];
        if (action_left > 2) action_left = action_left- 2;
        action_right = 3 - action_left; 
        
        if(pressed_left_prev == 1) {
          v_hybrid[action_left] =  v_hybrid[action_left] + st_par;
        } else {
           v_hybrid[action_right] =  v_hybrid[action_right] + st_par;
        }
       
      }
      
      level1_prob_choice2 = inv_logit( noise[i]*(v_hybrid[2]-v_hybrid[1]));
  
      log_lik[i] = log_lik[i] + bernoulli_lpmf( choice01 | level1_prob_choice2 );
      
       // generate posterior prediction for current trial
      y_pred[i,t] = bernoulli_rng(level1_prob_choice2);

      // log_lik[i] = log_lik[i] + categorical_logit_lpmf( choice | noise[i]*v_hybrid );
      // alternative model formulation
      // y_pred[i,t] = categorical_logit_rng( noise[i]*v_hybrid);

      // Observe Level2 and update Level1 of the chosen option
      //    dtQ(2) = subdata.points(t) - Q2(s2);                            % prediction error (2nd choice)
      delta = reward[i,t] - v_mf[state1[i,t], choice];
      if (delta >0) {
        a1_par = a1_g_par;
      } else {
        a1_par = a1_l_par;
      }
      
      v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + a1_par*delta;
        
      // After observing the reward at Level 2...
      // Update Level 2 v_mf of the chosen option. Level 2--> choose one of level 2 options and observe reward
      //    Q2(s2) = Q2(s2) + lr*dtQ(2);            
      v_2[state2[i,t]] =v_2[state2[i,t]] + a2_par*(reward[i,t] - v_2[state2[i,t]] );
      
      
       // forget others
      
      v_mf[state1[i,t], 3-choice] = (1-g_par)*v_mf[state1[i,t], 3-choice];
      v_mf[3-state1[i,t], choice] = (1-g_par)*v_mf[3-state1[i,t], choice];
      v_mf[3-state1[i,t], 3-choice] = (1-g_par )*v_mf[3-state1[i,t], 3-choice];
      
      // Update Level 1 v_mf with eligibility trace
       //    Qmf(s1,a) = Qmf(s1,a) + lambda*lr*dtQ(2);   
      // v_mf[state1[i,t], choice] = v_mf[state1[i,t], choice] + lambda[i] * a1[i] * (reward[i,t] - v_2[state2[i,t]] );
        } // end of t loop
        
    
        
    } // end of i loop
   } // end local
}

