Write ONE Python 3 file, no prose, no markdown fences: only code. UTF-8, `sys.stdout.reconfigure(encoding="utf-8")`. numpy allowed, scipy not.

TASK: port the Geant4 11.2.1 beta-minus spectrum of 40K (uniqueThirdForbidden) from the C++ excerpts below and compute
internal-bremsstrahlung (IB) yields using two beta shapes.

Constants: ME=511.0 keV; Q=1310.894 keV; daughter Z=20, A=40; CLHEP fine_structure_const = 1/137.03599976;
alphaZ = fsc*Z; Rnuc = 0.5*fsc*A**(1/3); V0 = 1.13*fsc**2*abs(Z)**(4/3); gamma0 = sqrt(1-alphaZ**2);
gc = [-0.1010678, 0.4245549, -0.6998588, 0.9512363, -0.5748646, 1.0].
Port EXACTLY (same formulas, same order, no simplification): Gamma(arg), ModSquared(re,im) with imMax=200, FermiFunction(W),
and ShapeFactor for case uniqueThirdForbidden (p_e = momentum in units of m_e, e_nu = (maxEnergy - ex) in units of m_e).
Spectrum on grid n=3000 midpoints of T in (0, Q): ex = T/ME; p = sqrt(ex*(ex+2)); en = Q/ME - ex;
f = p*(1+ex)*en*en * FermiFunction(1+ex) * ShapeFactor(p, en). Normalize so that np.trapezoid(f, T_keV) == 1. Call it (Tk, Pg).

Then `import ib_spectrum_kub as kub` (same folder, already exists; use kub.beta(E0, Z, shape) -> (T_keV, P normalized) and
kub.phi(W, k) where W = 1+T/ME and k in units of ME; do not rewrite them). Our shape: Tk2, Pk = kub.beta(Q, 20, "u3").
IB spectrum per beta decay for any (T, P): dNdk(k_keV) = np.trapezoid(P * kub.phi(1+T/ME, k_keV/ME) / ME, T).
ks = np.arange(0.5, 1310.0, 1.0) (keV). Print, in Russian, with sys.stdout UTF-8:
1) mean beta energy for Pg and Pk (trapezoid of T*P) and max relative difference of Pg vs np.interp(Tk, Tk2, Pk) in percent;
2) table for bands 25-150, 150-400, 400-1311, 25-1311 keV: sum of dNdk over ks in band for Geant4 shape, for our shape, and ratio;
3) save np.savetxt("../audit/gs45_k40_ib_shapes.csv" relative to the script file's folder, columns ks, dNdk_g4, dNdk_ours, header "k_keV,dNdk_g4shape,dNdk_ours_u3").
Put `if __name__ == "__main__": main()`. Guard division by zero in ratios.

C++ SOURCE (authoritative):
  // alphaZ = fine_structure_const*std::abs(Z);
  alphaZ = fine_structure_const*Z;

  // Nuclear radius in units of hbar/m_e/c
  G4double a13 = G4Pow::GetInstance()->Z13(A);
  Rnuc = 0.5*fine_structure_const*a13;

  // Electron screening potential in units of electron mass
  V0 = 1.13*fine_structure_const*fine_structure_const
           *std::pow(std::abs(Z), 4./3.);

  gamma0 = std::sqrt(1. - alphaZ*alphaZ);

  // Largest allowed value of im argument in ModSquared
//  imMax = std::log(DBL_MAX)/pi;
  imMax = 200.;   // actual value = 225.931, but use 200 to be safe
//  G4cout << " imMax = " << imMax << G4endl; 

  // Coefficients for gamma function with real argument
  gc[0] = -0.1010678;
  gc[1] =  0.4245549;
  gc[2] = -0.6998588;
  gc[3] =  0.9512363;
  gc[4] = -0.5748646;
  gc[5] = 1.0;
}


G4double G4BetaDecayCorrections::FermiFunction(const G4double& W)
{
  // Calculate the relativistic Fermi function.  Argument W is the
  // total electron energy in units of electron mass.

  G4double Wprime;
  if (Z < 0) {
    Wprime = W + V0;
  } else {
    Wprime = W - V0;
//    if (Wprime < 1.) Wprime = W;
    if (Wprime <= 1.00001) Wprime = 1.00001;
  }

  G4double p_e = std::sqrt(Wprime*Wprime - 1.);
  G4double eta = alphaZ*Wprime/p_e;
  G4double epieta = std::exp(pi*eta);
  G4double realGamma = Gamma(2.*gamma0+1);
  G4double mod2Gamma = ModSquared(gamma0, eta);

  // Fermi function
  G4double factor1 = 2*(1+gamma0)*mod2Gamma/realGamma/realGamma;
  G4double factor2 = epieta*std::pow(2*p_e*Rnuc, 2*(gamma0-1) );

  // Electron screening factor
  G4double factor3 = (Wprime/W)*std::sqrt( (Wprime*Wprime - 1.)/(W*W - 1.) );

  return factor1*factor2*factor3;
}


G4double
G4BetaDecayCorrections::ModSquared(const G4double& re, G4double im)
{
  // Calculate the squared modulus of the Gamma function 
  // with complex argument (re, im) using approximation B 
  // of Wilkinson, Nucl. Instr. & Meth. 82, 122 (1970).
  // Here, choose N = 1 in Wilkinson's notation for approximation B

  im = std::max(std::min(im, imMax), -imMax);
  G4double factor1 = std::pow( (1+re)*(1+re) + im*im, re+0.5);
  G4double factor2 = std::exp(2*im * std::atan(im/(1+re)));
  G4double factor3 = std::exp(2*(1+re));
  G4double factor4 = 2.*pi;
  G4double factor5 = std::exp( (1+re)/( (1+re)*(1+re) + im*im)/6 );
  G4double factor6 = re*re + im*im;
  return factor1*factor4*factor5/factor2/factor3/factor6;
}


G4double G4BetaDecayCorrections::Gamma(const G4double& arg)
{
  // Use recursion relation to get argument < 1
  G4double fac = 1.0;
  G4double x = arg - 1.;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl; 
  while (x > 1.0) { /* Loop checking, 01.09.2015, D.Wright */
    fac *= x;
    x -= 1.0;
    loop++;
    if (loop > 1000) {
      G4Exception("G4BetaDecayCorrections::Gamma()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
  }

  // Calculation of Gamma function with real argument
  // 0 < arg < 1 using polynomial from Abramowitz and Stegun
  G4double sum = gc[0];
  for (G4int i = 1; i < 6; i++) sum = sum*x + gc[i];

  return sum*fac;
}
    case (thirdForbidden) :
      break;

    case (uniqueThirdForbidden) :
      {
        G4double eta = alphaZ*std::sqrt(1. + p_e*p_e)/p_e;
        G4double gamma1 = std::sqrt(4. - alphaZ*alphaZ);
        G4double gamma2 = std::sqrt(9. - alphaZ*alphaZ);
        G4double gamma3 = std::sqrt(16. - alphaZ*alphaZ);
        G4double gamterm0 = Gamma(2.*gamma0+1.);
        G4double gamterm1 = gamterm0/Gamma(2.*gamma1+1.);
        G4double gamterm2 = gamterm0/Gamma(2.*gamma2+1.);
        G4double gamterm3 = gamterm0/Gamma(2.*gamma3+1.);

        G4double term1 = e_nu*e_nu*e_nu*e_nu*e_nu*e_nu*(1. + gamma0)/1260.;

        G4double term2 = 2.*(2. + gamma1)*e_nu*e_nu*e_nu*e_nu*p_e*p_e
                           *std::pow(twoPR, 2.*(gamma1-gamma0-1.) )
                           *gamterm1*gamterm1
                           *ModSquared(gamma1, eta)/ModSquared(gamma0, eta)/5.;

        G4double term3 = 60.*(3.+gamma2)*p_e*p_e*p_e*p_e*e_nu*e_nu
                             *std::pow(twoPR, 2.*(gamma2-gamma0-2.) )
                             *gamterm2*gamterm2
                             *ModSquared(gamma2, eta)/ModSquared(gamma0, eta);

        G4double term4 = 2240.*p_e*p_e*p_e*p_e*p_e*p_e*(4. + gamma3)
                             *std::pow(twoPR, 2.*(gamma3-gamma0-3.) )
                             *gamterm3*gamterm3
                             *ModSquared(gamma3, eta)/ModSquared(gamma0, eta);

