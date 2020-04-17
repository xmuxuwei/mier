#' @title Mie and Rayleigh function based on m_ri, lamda, diameter
#'
#' @description the funtions to calculate Q_ab, Q_ext, Q_sca,Q_pr, Q_back, g,Q_ratio
#'              the units of lamda and diameter should be nanometer,
#'              the Mie_Rayleigh function used a crossover of x == 0.01,
#'              when x < 0.01, the function returns results of Rayleigh scattering,
#'              when x > 0.01, the function returns Mie scattering.
#'              m_ri : refractive index,
#'              incident light wavelength: lamda
#'              particle diameter: diameter
#' @param
#'
#' @return
#'
#' @example mie(m_ri = 5+0.4i,lamda = 500, diameter = 300)
#'
#' @export


mie <- function(m_ri,lamda,diameter){
  #lamda the wavelength of indicent light, unit is nanometer
  #diameter, the diameter of the particles, unit is nanometer
  # m_ri is the complex refractive index of the particle,depending on chemical composition
  x = pi*diameter/lamda
  n_max = floor(x + 4*x^(1/3)+2)
  Q_sca= 0
  Q_back = 0
  Q_ext = 0
  g = 0
  i = 1
  while(i < n_max|i == n_max){
    Q_sca = Q_sca + 2/(x^2)*(2 *i +1)*(abs(an(m_ri,x,i))^2 +abs(bn(m_ri,x,i))^2)
    Q_back = Q_back + (2 * i + 1)*(-1)^(i) * (an(m_ri,x,i)- bn(m_ri,x,i))
    Q_ext = Q_ext + 2/(x^2)*(2 * i + 1)* Re(an(m_ri,x,i)+ bn(m_ri,x,i))

    asy1 =  (i * (i + 2)/(i+1))* Re( an(m_ri,x,n = i) * Conj(an(m_ri,x,n = i + 1)) +
                                       bn(m_ri,x,n = i) * Conj(bn(m_ri,x,n = i+1)))
    asy2 = (2 * i + 1)/(i * (i + 1))* Re(an(m_ri,x,n = i) * Conj(bn(m_ri,x,n = i)))
    g = g + (asy1 + asy2)
    i = i + 1}
  Q_back = (abs(Q_back))^2 / x^2
  Q_abs = Q_ext - Q_sca
  g = g*4/x^2/Q_sca
  Q_pr = Q_ext - g *Q_sca
  Q_ratio = Q_back / Q_sca

  return(tibble(Q_sca,Q_back,Q_ext,Q_abs,Q_pr,Q_ratio, g))}


Rayleigh <- function(m_ri,lamda,diameter){
  #computes the mie efficiencies of a spherical particles in the rayleigh regime( when x = pi*diameter/lamda << 1 )
  x = pi*diameter/lamda
  Q_sca = 8*x^4/3*(abs((m_ri^2 -1)/(m_ri^2+2)))^2
  Q_abs = 4 *x* Im((m_ri^2 -1)/(m_ri^2+2))
  Q_ext = Q_sca + Q_abs
  Q_back = 3 * Q_sca /2
  Q_ratio = 1.5
  Q_pr = Q_ext
  g = 0
  return(tibble(Q_sca,Q_back,Q_ext,Q_abs,Q_pr,Q_ratio,g))
}

Mie_Rayleigh <- function(m_ri,lamda,diameter,crossover = 0.01){
  #Mie-Rayleigh is auto function to choose mie or rayleigh based on crossover point 0.01
  # the crossover point is choose based on python package PyMieScatt (AutoMieQ) function
  x = pi*diameter/lamda
  if(x < crossover){
    return(Rayleigh(m_ri,lamda,diameter))
  } else{
    return(mie(m_ri,lamda,diameter))}
}
