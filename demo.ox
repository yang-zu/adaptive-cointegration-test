#include <oxstd.oxh>
#include "adaptcoint.ox"

main()
{
    decl n0 = 1000, dt = 1.0/n0;
    ranseed(20160119);

    decl rho = 0.4;
    decl sigbar = (1-rho)*unit(2) + rho*ones(2,2);
    // constant-volatility covariance path (3 x n0): [s11; s22; s12]
    decl s2 = (ones(1,n0)*sigbar[0][0]) | (ones(1,n0)*sigbar[1][1]) | (ones(1,n0)*sigbar[0][1]);
    decl s  = sig2tosig(s2);
    decl n = 500, step = n0/n, skip = range(0, n0-step, step);
    s = s[][skip];

    decl beta = <1;0>, psi = <0, 0.5 ; 0, 0.5>;

    decl Xnull = DGP(<0;0>,      beta, psi, s);   // r = 0 (no cointegration)
    decl Xalt  = DGP(<-15;0>/n,  beta, psi, s);   // cointegrated

    decl rn = coint_test(Xnull, 5);
    decl ra = coint_test(Xalt,  5);

    println("cols: k  trace  adaptLR  p_wb_trace  p_wb_adapt  p_vb_trace  p_vb_adapt");
    println("no cointegration : ", rn);
    println("cointegrated     : ", ra);
}
