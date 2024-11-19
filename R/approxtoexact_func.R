#' Convert the approximate allocation (proportion) to exact allocation (integer) without constraint
#'
#' @param n Sample size, must be a positive integer
#' @param w Approximate allocation/proportion, must be a real-valued vector, can get from running liftone_constrained_GLM or liftone_constrained_MLM
#'
#' @return allocation is the exact allocation or integer value of the number of subjects sampled from the group
#' @export
#'
#' @examples
#' exact_design = approxtoexact_func(n=600, w=c(0.2593526, 0.0000000, 0.0000000,
#' 0.1565024, 0.2891565, 0.0000000, 0.0000000, 0.2949885))
#'
#'

approxtoexact_func <- function(n, w) {
  # w[1:m] is a real-valued allocation, 0 <= w_i <= constraints[i]/n, sum_i w_i = 1
  # constraints[1:m] n_i <= constraints[i]
  m=length(w);
  constraints=rep(n,m);
  ans=floor(n*w);
  k=n-sum(ans);
  if(k>0) {
    stemp=(n*w-ans)*((ans+1) <= constraints);
    otemp=order(-stemp)[1:k];
    ans[otemp]=ans[otemp]+1;
  };
  #define S3 class
  output=list(allocation=ans);
  class(output)="list_output"
  return(output)
}
