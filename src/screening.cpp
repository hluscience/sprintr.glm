#include <cmath>
#include <queue>
#include <vector>
// we only include RcppArmadillo.h which automatically pulls Rcpp.h in for us
#include "RcppArmadillo.h"


// Define threshold and maximum iterations in Newton-Raphson
#define TH 1e-5
#define maxit 100000

// Logistic function
// [[Rcpp::export]]
arma::vec calc_mu_log(const arma::vec & x, const double & beta, const arma::vec & offset)
{
  return(1 / (1 + 1/exp(x*beta + offset)));
}

// Newton-Raphson for logistic
// [[Rcpp::export]]
double logistic(const arma::vec & x, const arma::vec & y, const arma::vec & offset){

  // setup
  double beta = 0.00;  // initial value for beta
  double diff = 10000; // initial value for diff
  int iter = 0; // iteration
  arma::vec mu;
  arma::vec W;
  double H;
  double U;
  double beta_change;
  double MML;

  // iterative
  while(diff > TH && iter < maxit){
    mu = calc_mu_log(x, beta, offset);
    W = mu%(1-mu);
    H = -arma::dot(x%W, x);
    U = arma::dot(x, (y - mu));

    // beta change
    beta_change = U/H;

    // update beta
    beta = beta - beta_change;

    // calculate how much we changed beta by in this iteration
    // if this is less than threshold, we'll break the while loop
    diff = beta_change*beta_change;

    iter = iter + 1;

    // if H is near to zero, beta_change approaches infinity,
    // so we need to break the while loop when H is small than TH
    if (abs(H) < TH){
      beta = 0;
      break;
    }

    if(iter == maxit) {
      printf("Convergence not reached after maxit = %d\n iterations in logistic screening step using newton-raphson for calculating MMLE", maxit);
    }
  }

  mu = calc_mu_log(x, beta, offset);

  // return marginal maximum likelihood
  MML = arma::dot(y, log(mu)) + arma::dot(1-y, log(1-mu));
  return MML;
}

// poisson function
// [[Rcpp::export]]
arma::vec calc_mu_pois(const arma::vec & x, const double & beta, const arma::vec & offset)
{
  return(exp(x*beta + offset));
}

// Newton-Raphson for poisson
// [[Rcpp::export]]
double poisson(const arma::vec & x, const arma::vec & y, const arma::vec & offset){

  // setup
  double beta = 0.00;  // initial value for beta
  double diff = 10000; // initial value for diff
  int iter = 0; // iteration
  arma::vec mu;
  double H;
  double U;
  double beta_change;
  double MML;

  // iterative
  while(diff > TH && iter < maxit){
    mu = calc_mu_pois(x,beta,offset);
    H = -arma::dot(x%mu, x);
    U = arma::dot(x, (y - mu));

    // beta change
    beta_change = U/H;

    // update beta
    beta = beta - beta_change;

    //calculate how much we changed beta by in this iteration
    //if this is less than threshold, we'll break the while loop
    diff = beta_change*beta_change;

    iter = iter + 1;

    // if H is near to zero, beta_change approaches infinity,
    // so we need to break the while loop when H is small than TH
    if (abs(H) < TH){
      beta = 0;
      break;
    }

    if(iter == maxit) {
      printf("Convergence not reached after maxit = %d\n iterations in poisson screening step using newton-raphson for calculating MMLE",maxit);
    }
  }

  mu = calc_mu_pois(x, beta, offset);
  MML = -sum(mu) + arma::dot(y, log(mu));
  return MML;
}

// multinomial function
// [[Rcpp::export]]
arma::mat calc_mu_multi(const arma::colvec & x, const arma::rowvec & beta, const arma::mat & offset)
{
  arma::mat exponent = exp(x*beta + offset);
  return(exponent.each_col()/arma::sum(exponent, 1));
}

// Newton-Raphson for multinomial
// [[Rcpp::export]]
double multi(const arma::vec & x, const arma::vec & y, const arma::mat & offset, const int & K){

  // setup
  int iter_outer = 0; // iteration for outer loop
  int iter_inner;
  double diff_outer = 10000; // initial value for diff for outer loop
  double diff_inner;
  arma::rowvec beta(K, arma::fill::zeros); // initial value for beta
  arma::rowvec beta_old;
  int k;
  double beta_k;
  arma::mat mu;
  arma::vec mu_k;
  arma::vec W_k;
  double H_k;
  double U_k;
  double beta_change_k;
  double MML = 0;

  while(diff_outer > TH && iter_outer < 100){
    beta_old = beta;
    for(k = 0; k < K-1; ++k){
      iter_inner = 0; // iteration for inner loop
      diff_inner = 10000; // initial value for diff for inner loop
      while(diff_inner > TH && iter_inner < maxit){
        // calculate probabilities using current estimate of beta
        beta_k = beta(k);
        mu =  calc_mu_multi(x, beta, offset); // Note: we can use the offset instead of offset-offset.col(K-1), because the offset.col(K-1) can be cancelled out by the exp term.
        mu_k = mu.col(k);

        // calculate W H U
        W_k = mu_k%(1-mu_k);
        H_k = -arma::dot(x%W_k, x);
        U_k = arma::dot(x, ((y == (k+1)) - mu_k));
        beta_change_k = U_k/H_k;

        // update beta
        beta_k = beta_k - beta_change_k;
        beta(k) = beta_k;

        // calculate how much we changed beta_k by in this iteration
        // if this is less than threshold, we'll break the while loop
        diff_inner = beta_change_k*beta_change_k;

        iter_inner = iter_inner + 1;

        // if H_k is near to zero, beta_change approaches infinity,
        // so we need to break the while loop when H_k is small than TH
        if (abs(H_k) < TH){
          beta_k = 0;
          beta(k) = beta_k;
          break;
        }

        if(iter_inner == maxit) {
          printf("Convergence not reached after maxit = %d\n inner iterations in multinomial screening step using newton-raphson Coordinate descent for calculating MMLE",maxit);
        }
      }
    }
    diff_outer = arma::max(abs(beta - beta_old));
    iter_outer = iter_outer + 1;
    if(iter_outer == maxit) {
      printf("Convergence not reached after maxit = %d\n outer iterations in multinomial screening step using newton-raphson Coordinate descent for calculating MMLE",maxit);
    }
  }

  // MML
  mu =  calc_mu_multi(x, beta, offset);
  for(k = 0; k < K; ++k){
    MML = MML + arma::dot((y == (k+1)), log(mu.col(k)));
  }

  return MML;
}

// delta function in ordinal
// [[Rcpp::export]]
arma::mat calc_delta_ordinal(const arma::mat & etaMat)
{
  return(1 / (1 + 1/exp(etaMat)));
}

// eta function in ordinal
// [[Rcpp::export]]
arma::rowvec calc_eta_ordinal(const arma::rowvec & delta)
{
  return(log(delta / (1 - delta)));
}

// span x matrix
// [[Rcpp::export]]
arma::vec tran_xspan(const arma::vec & x, const int & K){
  arma::vec xSpan = repelem(x, K-1, 1);
  return xSpan;
}

// span y matrix
// [[Rcpp::export]]
arma::mat tran_yspan(const arma::vec & y, const int & n, const int & K){

  arma::mat ySpan = arma::zeros(n, K);
  for (int i = 0; i < n; ++i){
    ySpan(i,y(i)-1) = 1;
  };
  return ySpan;
}


// Newton-Raphson for ordinal
// [[Rcpp::export]]
double ordinal(const arma::vec & x, const arma::vec & y, const arma::mat & offset, const int & n, const int & K){

  // setup
  double beta = 0.00;  // initial value for beta
  double diff = 10000; // initial value for diff
  int iter = 0; // iteration
  double score;
  double info;
  double beta_change;
  double MML;
  arma::vec xSpan;
  arma::mat ySpan;
  arma::mat ttipDiag;
  arma::mat ttipOffDiag;
  arma::mat ttip;
  arma::mat etaMat;
  arma::mat deltaMat;
  arma::mat pMat;
  arma::mat pkplusone;
  arma::mat pMatFull;
  arma::vec d;
  arma::mat uMat;
  arma::mat uminusdMat;
  arma::mat wpMat;
  arma::vec wpkplusone;
  arma::vec xEach;
  arma::rowvec uminusd;
  arma::rowvec delta;
  arma::mat q;
  arma::mat tt;
  arma::mat sigInv;
  arma::mat w;

  // span x and y matrix
  xSpan = tran_xspan(x, K);
  ySpan = tran_yspan(y, n, K);

  ttipDiag.eye(K-1, K-1);
  ttipOffDiag = join_rows(-ttipDiag.cols(arma::span(1,K-2)), arma::zeros(K-1,1));
  ttip = ttipDiag + ttipOffDiag;

  while(diff > TH && iter < maxit){
    // etaMat -> deltaMat -> pMat -> pMatFull
    etaMat = arma::trans(arma::reshape(xSpan*beta, K-1, n)) + offset; // linear part
    deltaMat = calc_delta_ordinal(etaMat); // use inverse link
    pMat = deltaMat - join_rows(arma::zeros(n,1), deltaMat.cols(arma::span(0,K-3)));
    pkplusone = arma::ones(n,1) - arma::sum(pMat, 1);
    pMatFull = join_rows(pMat, pkplusone);

    // prepare for calcualting score and info
    d = ySpan.col(K-1) / pkplusone;
    uMat = ySpan.cols(arma::span(0,K-2)) / pMat;
    uminusdMat = uMat - arma::repmat(d, 1, K-1);

    wpMat = arma::ones(n, K-1) / pMat;
    wpkplusone = arma::ones(n, 1) / pkplusone;

    score = 0.00;
    info = 0.00;
    for (int i = 0; i < n; ++i){
      // compute score term
      xEach = xSpan.rows(arma::span((K-1)*i, (K-1)*i+K-2));
      uminusd = uminusdMat.row(i);
      delta = arma::cumsum(pMat.row(i));
      q = reshape(vectorise(repelem(delta%(1-delta), K-1, 1))%vectorise(ttip), K-1, K-1);
      score = score + as_scalar(arma::trans(xEach) * (arma::trans(q) * arma::trans(uminusd)));

      // compute info term
      sigInv = arma::diagmat(wpMat.row(i)) + wpkplusone(i);
      w = arma::trans(q) * sigInv * q;
      info = info + as_scalar(arma::trans(xEach) *(w * xEach));
    };

    // beta change
    beta_change = score/info;

    // update beta
    beta = beta + beta_change;

    //calculate how much we changed beta by in this iteration
    //if this is less than threshold, we'll break the while loop
    diff = beta_change*beta_change;

    iter = iter + 1;

    if(iter == maxit) {
      printf("Convergence not reached after maxit = %d\n iterations in ordinal screening step using newton-raphson for calculating MMLE",maxit);
    }
  }

  // etaMat -> deltaMat -> pMat -> pMatFull -> MML
  etaMat = arma::trans(arma::reshape(xSpan*beta, K-1, n)) + offset;
  deltaMat = calc_delta_ordinal(etaMat);
  pMat = deltaMat - join_rows(arma::zeros(n,1), deltaMat.cols(arma::span(0,K-3)));
  pkplusone = arma::ones(n,1) - arma::sum(pMat, 1);
  pMatFull = join_rows(pMat, pkplusone);
  MML = arma::accu(ySpan % log(pMatFull));

  return MML;
}


// Newton-Raphson for ordinal sparse
// [[Rcpp::export]]
double ordinal_sparse(const arma::vec x, const arma::vec & y, const arma::sp_mat & yMat, const arma::mat & offset, const int & n, const int & K){

  // extract subvec and submat with non-zero elements
  arma::uvec idx = arma::find(x != 0);
  const int n_sub = idx.size();
  arma::vec x_sub = x.elem(idx);
  arma::vec y_sub = y.elem(idx);
  arma::mat offset_sub = offset.rows(idx);

  // setup
  double beta = 0.00;  // initial value for beta
  double diff = 10000; // initial value for diff
  int iter = 0; // iteration
  double score;
  double info;
  double beta_change;
  double MML;
  arma::vec xSpan;
  arma::mat ySpan;
  arma::mat ttipDiag;
  arma::mat ttipOffDiag;
  arma::mat ttip;
  arma::mat etaMat;
  arma::mat deltaMat;
  arma::mat pMat;
  arma::mat pkplusone;
  arma::mat pMatFull;
  arma::vec d;
  arma::mat uMat;
  arma::mat uminusdMat;
  arma::mat wpMat;
  arma::vec wpkplusone;
  arma::rowvec uminusd;
  arma::mat q;
  arma::mat sigInv;
  arma::mat w;

  // span x and y matrix
  xSpan = tran_xspan(x_sub, K);
  ySpan = tran_yspan(y_sub, n_sub, K);

  ttipDiag.eye(K-1, K-1);
  ttipOffDiag = join_rows(-ttipDiag.cols(arma::span(1,K-2)), arma::zeros(K-1,1));
  ttip = ttipDiag + ttipOffDiag;

  while(diff > TH && iter < maxit){
    // etaMat -> deltaMat -> pMat -> pMatFull
    etaMat = arma::trans(arma::reshape(xSpan*beta, K-1, n_sub)) + offset_sub;
    deltaMat = calc_delta_ordinal(etaMat);
    pMat = deltaMat - join_rows(arma::zeros(n_sub,1), deltaMat.cols(arma::span(0,K-3)));
    pkplusone = arma::ones(n_sub,1) - arma::sum(pMat, 1);
    pMatFull = join_rows(pMat, pkplusone);

    // prepare for calcualting score and info
    d = ySpan.col(K-1) / pkplusone;
    uMat = ySpan.cols(arma::span(0,K-2)) / pMat;
    uminusdMat = uMat - arma::repmat(d, 1, K-1);

    wpMat = arma::ones(n_sub, K-1) / pMat;
    wpkplusone = arma::ones(n_sub, 1) / pkplusone;

    score = 0.00;
    info = 0.00;
    for (int i = 0; i < n_sub; ++i){
      // compute score term
      uminusd = uminusdMat.row(i);
      q = reshape(vectorise(repelem(deltaMat.row(i)%(1-deltaMat.row(i)), K-1, 1))%vectorise(ttip), K-1, K-1);
      score = score + x_sub(i) * arma::accu(arma::trans(q) * arma::trans(uminusd));
      // compute info term
      sigInv = arma::diagmat(wpMat.row(i)) + wpkplusone(i);
      w = arma::trans(q) * sigInv * q;
      info = info + arma::accu(w) * x_sub(i) * x_sub(i);
    };

    // if info is near to zero, beta_change approaches infinity,
    // so we need to break the while loop when info is small than TH
    if (abs(info) < TH){
      beta = 0.00;
      break;
    }

    // beta change
    beta_change = score/info;

    // update beta
    beta = beta + beta_change;

    //calculate how much we changed beta by in this iteration
    //if this is less than threshold, we'll break the while loop
    diff = beta_change*beta_change;

    iter = iter + 1;

    if(iter == maxit) {
      printf("Convergence not reached after maxit = %d\n iterations in screening step using poisson newton-raphson for calculating MMLE",maxit);
    }
  }

  // etaMat -> deltaMat -> pMat -> pMatFull -> MML
  etaMat = arma::trans(arma::reshape(tran_xspan(x, K)*beta, K-1, n)) + offset;
  deltaMat = calc_delta_ordinal(etaMat);
  pMat = deltaMat - join_rows(arma::zeros(n,1), deltaMat.cols(arma::span(0,K-3)));
  pkplusone = arma::ones(n,1) - arma::sum(pMat, 1);
  pMatFull = join_rows(pMat, pkplusone);
  MML = arma::accu(yMat % log(pMatFull));

  return MML;
}

class Triplet
{
private:
  // j, k are the indices of the interaction
  int j, k;
  // the MML of the j,k interaction with some target
  double MML;
public:
  Triplet(int _j, int _k, double _MML)
  {
    j = _j;
    k = _k;
    MML = _MML;
  }
  int get_j() const{
    return j;
  }
  int get_k() const{
    return k;
  }
  double get_MML() const{
    return MML;
  }
};

// user specified comparator for MML
class Triplet_Comparator{
public:
  bool operator() (const Triplet & p1, const Triplet & p2) const {
    // note that here we use > instead of <
    // so that the resulting heap is a min-heap
    // (by default the C++ priority_queue is a max-heap)
    return p1.get_MML() > p2.get_MML();
  }
};

//' Sure Independence Screening in Step 2
//'
//' @param x a n-by-p matrix of main effects, with i.i.d rows, and each row represents a vector of observations of p main-effects
//' @param y a response vector of size \code{n}. \code{y} is a vector of continuous data in "gaussian" family; a vector of positive count data in "poisson" family; a vector of 0-1 binary data in "binomial" family. Notice that \code{y} should be a factor when using the "ordinal" family for ordinal response.
//' @param nlev the number of levels for response variable. Set to 1 for "gaussian" and "poisson" family.
//' @param offset the linear predictions from Step 1.
//' @param num_keep the number of candidate interactions in Step 2. Default to be n.
//' @param square An indicator of whether squared effects should be considered in Step 1 (NOT Step 2!). square == TRUE if squared effects have been considered in Step 1, i.e., squared effects will NOT be considered in Step 2.
//' @param main_effect An indicator of whether main effects should also be screened. Default to be false.
//' @param family Either a character string representing one of the built-in exponential families, including "gaussian", "poisson", "binomial", "ordinal". Default is "gaussian".
//' @return a three-column matrix representing the index pair of the selected interactions and its MML with the response vector given by the \code{offset}.
//' @export
// [[Rcpp::export]]
arma::mat screen_cpp(const arma::mat & x, const arma::vec & y, const int & nlev, const arma::mat & offset, const int & num_keep, const bool square = false, const bool main_effect = false, std::string family = "gaussian"){
 // main is TRUE if we want to screen main effects simultaneously
 // square is TRUE if squared effects have been considered in the first step
 // i.e., squared effects will NOT be considered in the screening step!
 // screening all the interaction x_j*x_k based on its MML
 // return the indices (j, k) of the top k interactions

 const int n = x.n_rows;
 const int p = x.n_cols;
 // j, k are used for indexing interactions
 int j, k;
 int count = 0;
 double cur_MML = 0.0;
 arma::vec temp_vec(n);
 // matrix of num_keep by 2, storing the results
 arma::mat result = arma::zeros<arma::mat>(num_keep, 3);

 // Creates a Min heap of Triplets (ordered by MML)
 std::priority_queue <Triplet, std::vector<Triplet>, Triplet_Comparator> min_heap;

 // squared effects are considered in the first step
 // not include them in the screening
 if(square){
   for(j = 0; j < p-1; ++j){
     for (k = j+1; k < p; ++k){
       temp_vec = x.col(j) % x.col(k);
       temp_vec = (temp_vec - arma::mean(temp_vec))/arma::stddev(temp_vec);
       if(family == "gaussian"){
         // correlation is equivalent to MMLE in gaussian case
         // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
         cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
       }else if (family == "binomial"){
         cur_MML = logistic(temp_vec, y, offset);
       }else if(family == "poisson"){
         cur_MML = poisson(temp_vec, y, offset);
       }else if(family == "multinomial"){
         cur_MML = multi(temp_vec, y, offset, nlev);
       }else if(family == "ordinal"){
         cur_MML = ordinal(temp_vec, y, offset, n, nlev);
       }
       if (count < num_keep){
         // use emplace instead of push
         // emplace avoids constructing a temporary Triplet
         // emplace requires C++11, to use it, add
         // // [[Rcpp::plugins(cpp11)]] in the top of the file
         // min_heap.emplace(j, k, cur_MML);
         min_heap.push(Triplet(j, k, cur_MML));
         ++count;
       } else{
         if (cur_MML > min_heap.top().get_MML()){
           // pop the current top element
           min_heap.pop();
           // insert the (j,k) interaction
           //min_heap.emplace(j, k, cur_MML);
           min_heap.push(Triplet(j, k, cur_MML));
         }
       }
     }
   }
 }
 else{
   // squared effects are not considered in the first step
   // include them in the screening
   for(j = 0; j < p; ++j){
     for (k = j; k < p; ++k){
       temp_vec = x.col(j) % x.col(k);
       temp_vec = (temp_vec - arma::mean(temp_vec))/arma::stddev(temp_vec);
       if(family == "gaussian"){
         // correlation is equivalent to MMLE in gaussian case
         // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
         cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
       }else if(family == "binomial"){
         cur_MML = logistic(temp_vec, y, offset);
       }else if(family == "poisson"){
         cur_MML = poisson(temp_vec, y, offset);
       }else if(family == "multinomial"){
         cur_MML = multi(temp_vec, y, offset, nlev);
       }else if(family == "ordinal"){
         cur_MML = ordinal(temp_vec, y, offset, n, nlev);
       }
       if (count < num_keep){
         // use emplace instead of push
         // emplace avoids constructing a temporary Triplet
         // emplace requires C++11, to use it, add
         // // [[Rcpp::plugins(cpp11)]] in the top of the file
         // min_heap.emplace(j, k, cur_MML);
         min_heap.push(Triplet(j, k, cur_MML));
         ++count;
       } else{
         if (cur_MML > min_heap.top().get_MML()){
           // pop the current top element
           min_heap.pop();
           // insert the (j,k) interaction
           //min_heap.emplace(j, k, cur_MML);
           min_heap.push(Triplet(j, k, cur_MML));
         }
       }
     }
   }
 }

 // if we need to screen main effects also
 if(main_effect){
   for(k = 0; k < p; ++k){
     temp_vec = x.col(k);
     temp_vec = (temp_vec - arma::mean(temp_vec))/arma::stddev(temp_vec);
     if(family == "gaussian"){
       // correlation is equivalent to MMLE in gaussian case
       // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
       cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
     }else if(family == "binomial"){
       cur_MML = logistic(temp_vec, y, arma::zeros(n));
     }else if(family == "poisson"){
       cur_MML = poisson(temp_vec, y, arma::zeros(n));
     }else if(family == "multinomial"){
       cur_MML = multi(temp_vec, y, arma::zeros(n, nlev), nlev);
     }else if(family == "ordinal"){
       cur_MML = ordinal(temp_vec, y, arma::zeros(n, nlev), n, nlev);
     }
     if (cur_MML > min_heap.top().get_MML()){
       // pop the current top element
       min_heap.pop();
       // insert the (j,k) interaction
       // for main effects, it corresponds to (0, k)
       min_heap.push(Triplet(-1, k, cur_MML));
     }
   }
 }

 // now store the selected indices of interactions
 for (count = 0; count < num_keep; ++count){
   result(count, 0) = min_heap.top().get_j() + 1;
   result(count, 1) = min_heap.top().get_k() + 1;
   result(count, 2) = min_heap.top().get_MML();
   min_heap.pop();
 }
 return result;
}

//' Sure Independence Screening in Step 2 for sparse design matrix
//'
//' @param x a n-by-p sparse matrix of main effects.
//' @param y a response vector of size \code{n}. \code{y} is a vector of continuous data in "gaussian" family; a vector of positive count data in "poisson" family; a vector of 0-1 binary data in "binomial" family. Notice that \code{y} should be a factor when using the "ordinal" family for ordinal response.
//' @param yMat Default to be empty in "gaussian", "binomial" and "poisson" family. In "ordinal" family, it transforms an ordinal response vector \code{y} to a sparse \code{n}-by-\code{nlev} matrix with element 1 denoting the i-th observation with the j-th ordered level.
//' @param nlev the number of levels for response variable. Set to 1 for "gaussian" and "poisson" family.
//' @param offset the linear predictions from Step 1.
//' @param num_keep the number of candidate interactions in Step 2. Default to be n.
//' @param square An indicator of whether squared effects should be considered in Step 1 (NOT Step 2!). square == TRUE if squared effects have been considered in Step 1, i.e., squared effects will NOT be considered in Step 2.
//' @param main_effect An indicator of whether main effects should also be screened. Default to be false.
//' @param family Either a character string representing one of the built-in exponential families, including "gaussian", "poisson", "binomial", "ordinal". Default is "gaussian".
//' @return a three-column matrix representing the index pair of the selected interactions and its MML with the response vector given by the \code{offset}.
//' @export
// [[Rcpp::export]]
arma::mat screen_sparse_cpp(const arma::sp_mat & x, const arma::vec & y, const arma::sp_mat & yMat, const int & nlev, const arma::mat & offset, const int & num_keep, const bool square = false, const bool main_effect = false, std::string family = "gaussian"){
 // main is TRUE if we want to screen main effects simultaneously
 // square is TRUE if squared effects have been considered in the first step
 // i.e., squared effects will NOT be considered in the screening step!
 // screening all the interaction x_j*x_k based on its MML
 // return the indices (j, k) of the top k interactions

 const int n = x.n_rows;
 const int p = x.n_cols;
 // j, k are used for indexing interactions
 int j, k;
 int count = 0;
 double cur_MML = 0.0;
 arma::vec temp_vec(n);
 // matrix of num_keep by 2, storing the results
 arma::mat result = arma::zeros<arma::mat>(num_keep, 3);

 // Creates a Min heap of Triplets (ordered by MML)
 std::priority_queue <Triplet, std::vector<Triplet>, Triplet_Comparator> min_heap;

 // squared effects are considered in the first step
 // not include them in the screening
 if(square){
   for(j = 0; j < p-1; ++j){
     for (k = j+1; k < p; ++k){
       temp_vec = x.col(j) % x.col(k);
       if(family == "gaussian"){
         // correlation is equivalent to MMLE in gaussian case
         // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
         cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
       }else if(family == "binomial"){
         cur_MML = logistic(temp_vec, y, offset);
       }else if(family == "poisson"){
         cur_MML = poisson(temp_vec, y, offset);
       }else if(family == "multinomial"){
         cur_MML = multi(temp_vec, y, offset, nlev);
       }else if(family == "ordinal"){
         cur_MML = ordinal_sparse(temp_vec, y, yMat, offset, n, nlev);
       }
       if (count < num_keep){
         // use emplace instead of push
         // emplace avoids constructing a temporary Triplet
         // emplace requires C++11, to use it, add
         // // [[Rcpp::plugins(cpp11)]] in the top of the file
         // min_heap.emplace(j, k, cur_MML);
         min_heap.push(Triplet(j, k, cur_MML));
         ++count;
       } else{
         if (cur_MML > min_heap.top().get_MML()){
           // pop the current top element
           min_heap.pop();
           // insert the (j,k) interaction
           //min_heap.emplace(j, k, cur_MML);
           min_heap.push(Triplet(j, k, cur_MML));
         }
       }
     }
   }
 }
 else{
   // squared effects are not considered in the first step
   // include them in the screening
   for(j = 0; j < p; ++j){
     for (k = j; k < p; ++k){
       temp_vec = x.col(j) % x.col(k);
       if(family == "gaussian"){
         // correlation is equivalent to MMLE in gaussian case
         // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
         cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
       }else if(family == "binomial"){
         cur_MML = logistic(temp_vec, y, offset);
       }else if(family == "poisson"){
         cur_MML = poisson(temp_vec, y, offset);
       }else if(family == "multinomial"){
         cur_MML = multi(temp_vec, y, offset, nlev);
       }else if(family == "ordinal"){
         cur_MML = ordinal_sparse(temp_vec, y, yMat, offset, n, nlev);
       }
       if (count < num_keep){
         // use emplace instead of push
         // emplace avoids constructing a temporary Triplet
         // emplace requires C++11, to use it, add
         // // [[Rcpp::plugins(cpp11)]] in the top of the file
         // min_heap.emplace(j, k, cur_MML);
         min_heap.push(Triplet(j, k, cur_MML));
         ++count;
       } else{
         if (cur_MML > min_heap.top().get_MML()){
           // pop the current top element
           min_heap.pop();
           // insert the (j,k) interaction
           //min_heap.emplace(j, k, cur_MML);
           min_heap.push(Triplet(j, k, cur_MML));
         }
       }
     }
   }
 }

 // if we need to screen main effects also
 if(main_effect){
   for(k = 0; k < p; ++k){
     temp_vec = x.col(k);
     if(family == "gaussian"){
       // correlation is equivalent to MMLE in gaussian case
       // we use cur_MML instead of cur_MMLE to simplify notation in min_heap
       cur_MML = fabs(arma::as_scalar(arma::cor(temp_vec, y)));
     }else if(family == "binomial"){
       cur_MML = logistic(temp_vec, y, arma::zeros(n));
     }else if(family == "poisson"){
       cur_MML = poisson(temp_vec, y, arma::zeros(n));
     }else if(family == "multinomial"){
       cur_MML = multi(temp_vec, y, arma::zeros(n, nlev), nlev);
     }else if(family == "ordinal"){
       cur_MML = ordinal_sparse(temp_vec, y, yMat, arma::zeros(n, nlev), n, nlev);
     }

     if (cur_MML > min_heap.top().get_MML()){
       // pop the current top element
       min_heap.pop();
       // insert the (j,k) interaction
       // for main effects, it corresponds to (0, k)
       min_heap.push(Triplet(-1, k, cur_MML));
     }
   }
 }

 // now store the selected indices of interactions
 for (count = 0; count < num_keep; ++count){
   result(count, 0) = min_heap.top().get_j() + 1;
   result(count, 1) = min_heap.top().get_k() + 1;
   result(count, 2) = min_heap.top().get_MML();
   min_heap.pop();
 }
 return result;
}
