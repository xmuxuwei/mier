#' @title the an,bn,cn and dn function
#'
#' @description The functions were used to calculate internal and external field coefficients,
#'               m_ri is the complex refractive index,
#'                x is size parameter, n is order parameter,
#'               these functions are the basics for mie calculations.
#' @param
#'
#' @return complex
#'
#' @example an(m_ri = 5+0.4i)
#'
#' @export

j <- function(n,z) {sqrt(0.5*pi/z) * BesselJ(z,  nu = n + 0.5)}

y <- function(n,z) {sqrt(0.5*pi/z) * BesselY(z,  nu = n + 0.5)}

h1 <- function(m, n, z) {sqrt(0.5*pi/z) * BesselH(m, z, nu = n + 0.5)}

an <- function(m_ri, x, n) {
  (m_ri^2 * j(z = m_ri * x, n) * (x * j(z = x, n-1) - n*j(z = x,n)) -
     j(z = x , n) * (m_ri * x * j(z = m_ri*x,  n-1) - n*j(z = m_ri* x,  n)))/
    (m_ri^2 * j(z = m_ri * x,  n) * (x * h1(m = 1,z = x, n-1)- n * h1(m = 1,z = x,  n)) -
       h1(m = 1,z = x, n) * (m_ri*x * j(z = m_ri*x,  n-1) - n * j(z = m_ri* x,  n)))}

bn <- function(m_ri, x, n) {
  (j(z = m_ri*x,  n) * (x * j(z = x,  n-1) - n * j(z =x, n)) -
     j(z = x, n) * (m_ri*x * j(z = m_ri * x,  n-1) - n * j(z = m_ri* x, n)))/
    (j(z = m_ri*x,  n) * (x* h1(m = 1,z = x, n-1)- n * h1(m=1,z = x, n))-
       h1(m = 1,z = x, n) * (m_ri*x * j(z = m_ri*x,  n-1) - n * j(z = m_ri* x, n)))}


cn <- function(m_ri, x, n) {
  (j(z = x,  n)* (x * h1(m=1,z= x, n-1)-n*h1(m = 1,z = x,  n))-
     h1(m = 1,z = x, n) * (x * j(z = x,  n-1) - n*j(z = x,  n)))/
    (j(z = m_ri*x,  n) * (x*h1(m =1,z = x, n-1)-n*h1(m = 1,z = x,  n))-
       h1(m = 1,z = x, n) * (m_ri * x * j(z = m_ri*x,  n-1) - n*j(z = m_ri* x,  n)))}


dn <- function(m_ri, x, n) {
  (m_ri* j(z =x,n) *(x*h1(m=1,z= x, n-1)-n*h1(m=1,z = x, n)) -
     m_ri* h1(m=1,z =x,n) *(x * j(z=x,  n-1) - n*j(z =x, n)))/
    (m_ri^2 * j(z = m_ri * x,  n) * (x * h1(m=1,z= x, n-1)-n* h1(m = 1,z = x, n)) -
       h1(m = 1,z = x, n) * (m_ri * x * j(z= m_ri*x,  n-1) - n* j(z = m_ri* x, n)))}
