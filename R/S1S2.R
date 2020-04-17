#' @title S1S2
#'
#' @description the function calculates the scattering amplitudes S1 and S2,
#'              base on input angle:theta,
#'              refractive index : m_ri,
#'              incident light wavelength:lamda,
#'              partile diameter: diameter
#'
#'
#' @param m_ri, x, n
#'
#' @return
#'
#' @example mie(m_ri = 5+0.4i,lamda = 500, diameter = 300)
#'
#' @export
#'
pi_n <- function(n,theta){
  if(n ==0){return(0)}else{
    if(n == 1) {
      return(1)
    } else {
      i = 2
      pi_1 = 1
      pi_0 = 0
      while (i < n | i == n) {
        pi = (2*i-1)/(i-1)*cos(theta)*pi_1 - i/(i-1)*pi_0
        pi_0 = pi_1
        pi_1 = pi
        i= i+1
      }
      return(pi)
    }
  }
}

tau_n <-function(n,theta){
  if(n==0){return(0)}
  else{n * cos(theta)* pi_n(n=n,theta) - (n+1) * pi_n(n-1,theta)}}


mie_S1S2 <- function(theta,m_ri,lamda,diameter){
  #S1 and S2 are the scattering amplitudes as a function of theta, m_ri and x
  x = 2* pi* diameter / lamda
  n_max = floor(x + 4*x^(1/3) + 2)
  i = 1
  S1 = 0
  S2 = 0

  while(i < n_max |i == n_max){
    S1 = S1 + {(2*n+1)/(n*(n+1))*
        (an(m_ri,x,n)* pi_n(n, theta) + bn(m_ri,x,n)*tau_n(n,theta))}
    S2 = S2 +  {(2*n+1)/(n*(n+1))*
        (an(m_ri,x,n)*tau_n(n,theta) + bn(m_ri,x,n)* pi_n(n, theta))}
    i = i+1}
  return(tibble(S1,S2))
}
