
// Discontinuity sensor
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Real disconSensor(Real pp, Real pl,
                                                      Real pr) {
  Real pjst = pr + 2.0_rt * pp + pl;
  Real ptvd = std::abs(pr - pp) + std::abs(pp - pl);
  return std::abs(2.0_rt * (pr - 2.0_rt * pp + pl) /
                  (pjst + ptvd + Real(1.0e-40)));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ComputeSensorLambda(
    int i, int j, int k, const auto& prims, const auto& lambda, const auto& sen,
    const PROB::ProbClosures& closures) {
  // Real pp=prims(i,j,k,QPRES)*prims(i,j,k,QRHO);
  // sen(i,j,k,0) =
  // disconSensor(pp,prims(i-1,j,k,QPRES)*prims(i-1,j,k,QRHO),prims(i+1,j,k,QPRES)*prims(i+1,j,k,QRHO));
  // sen(i,j,k,1) =
  // disconSensor(pp,prims(i,j-1,k,QPRES)*prims(i,j-1,k,QRHO),prims(i,j+1,k,QPRES)*prims(i,j+1,k,QRHO));
  // sen(i,j,k,2) =
  // disconSensor(pp,prims(i,j,k-1,QPRES)*prims(i,j,k-1,QRHO),prims(i,j,k+1,QPRES)*prims(i,j,k+1,QRHO));

  // Real pp=prims(i,j,k,QRHO);
  // sen(i,j,k,0) = disconSensor(pp,prims(i-1,j,k,QRHO),prims(i+1,j,k,QRHO));
  // sen(i,j,k,1) = disconSensor(pp,prims(i,j-1,k,QRHO),prims(i,j+1,k,QRHO));
  // sen(i,j,k,2) = disconSensor(pp,prims(i,j,k-1,QRHO),prims(i,j,k+1,QRHO));

  Real pp = prims(i, j, k, QPRES);
  sen(i, j, k, 0) =
      disconSensor(pp, prims(i - 1, j, k, QPRES), prims(i + 1, j, k, QPRES));
  sen(i, j, k, 1) =
      disconSensor(pp, prims(i, j - 1, k, QPRES), prims(i, j + 1, k, QPRES));
  sen(i, j, k, 2) =
      disconSensor(pp, prims(i, j, k - 1, QPRES), prims(i, j, k + 1, QPRES));

  Real ux = prims(i, j, k, QU);
  Real uy = prims(i, j, k, QV);
  Real uz = prims(i, j, k, QW);
  Real cs = sqrt(closures.gamma * prims(i, j, k, QPRES) / prims(i, j, k, QRHO));
  lambda(i, j, k, 0) =
      std::abs(ux) + cs;  // max(std::abs(ux+cs),std::abs(ux-cs));
  lambda(i, j, k, 1) =
      std::abs(uy) + cs;  // max(std::abs(uy+cs),std::abs(uy-cs));
  lambda(i, j, k, 2) =
      std::abs(uz) + cs;  // max(std::abs(uz+cs),std::abs(uz-cs));
  // lambda(i,j,k,0) = max(std::abs(ux+cs),std::abs(ux-cs));
  // lambda(i,j,k,1) = max(std::abs(uy+cs),std::abs(uy-cs));
  // lambda(i,j,k,2) = max(std::abs(uz+cs),std::abs(uz-cs));
}

// calculates dissipative flux at i-1/2, j-1/2 and k-1/2
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void JSTflux(
    int i, int j, int k, int n, const auto& lambda, const auto& sensor,
    const auto& cons, const auto& nfabfx, const auto& nfabfy,
    const auto& nfabfz, const PROB::ProbClosures& closures) {
  Real dw, sen, rr, fdamp;
  Real u_ijk = cons(i, j, k, n);

  // x-dir
  dw = (u_ijk - cons(i - 1, j, k, n));
  sen = closures.Cshock * max(sensor(i - 1, j, k, 0), sensor(i, j, k, 0));
  rr = max(lambda(i - 1, j, k, 0), lambda(i, j, k, 0));
  fdamp = cons(i + 1, j, k, n) - 3.0 * cons(i, j, k, n) +
          3.0 * cons(i - 1, j, k, n) - cons(i - 2, j, k, n);
  nfabfx(i, j, k, n) -=
      (sen * dw - max(0.0, closures.Cdamp - sen) * fdamp) * rr;

  // y-dir
  dw = (u_ijk - cons(i, j - 1, k, n));
  sen = closures.Cshock * max(sensor(i, j - 1, k, 0), sensor(i, j, k, 0));
  rr = max(lambda(i, j - 1, k, 0), lambda(i, j, k, 0));
  fdamp = cons(i, j + 1, k, n) - 3.0 * cons(i, j, k, n) +
          3.0 * cons(i, j - 1, k, n) - cons(i, j - 2, k, n);
  nfabfy(i, j, k, n) -=
      (sen * dw - max(0.0, closures.Cdamp - sen) * fdamp) * rr;

  // z-dir
  dw = (u_ijk - cons(i, j, k - 1, n));
  sen = closures.Cshock * max(sensor(i, j, k - 1, 0), sensor(i, j, k, 0));
  rr = max(lambda(i, j, k - 1, 0), lambda(i, j, k, 0));
  fdamp = cons(i, j, k - 1, n) - 3.0 * cons(i, j, k, n) +
          3.0 * cons(i, j, k - 1, n) - cons(i, j, k - 1, n);
  nfabfz(i, j, k, n) -=
      (sen * dw - max(0.0, closures.Cdamp - sen) * fdamp) * rr;
}

