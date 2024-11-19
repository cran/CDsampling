#' Find (constrained) uniform exact allocation of the study for bounded design
#'
#' @param Ni a vector with size m, upper bound for exact design of each category/stratification group, if unconstrained, use Inf vector, the function will return unbounded uniform allocation
#' @param nsample a number, the sample size
#'
#' @return n is the constrained/unconstrained uniform exact allocation
#' @export
#'
#' @examples
#'
#' bounded_uniform(Ni=c(50, 40, 10, 200, 150, 50), nsample=200)
#'

bounded_uniform = function(Ni, nsample){
  N = sum(Ni) #population size
  nsample.temp = nsample #sample size
  m = m.temp = length(Ni) #num of categories
  Ni.temp = Ni #used in loop
  n = rep(NA, m) #allocation

  while(floor(nsample.temp/m.temp)>min(Ni.temp)){
    i = which.min(Ni.temp)
    n[i]=Ni[i]
    nsample.temp = nsample.temp - Ni[i]
    m.temp = m.temp-1
    Ni.temp[i] = Inf
  }

  k = floor(nsample.temp/m.temp)

  n = ifelse(is.na(n), k, n)

  diff = nsample - sum(n)
  if(diff>0){
    id = sample(seq(1,m), diff)
    n[id]=n[id]+1;
  };
  #define S3 class
  output <- list(allocation=n)
  class(output)<-"list_output"
  return(output)

}
