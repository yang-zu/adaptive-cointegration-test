# Adaptive cointegration test with nonstationary volatility — Ox

Ox implementation of the adaptive cointegration test proposed in:

> Boswijk, H. P. and Zu, Y. (2022). Adaptive testing for cointegration with
> nonstationary volatility. *Journal of Business and Economic Statistics*, 40, 744–755.
> [doi:10.1080/07350015.2020.1867558](https://doi.org/10.1080/07350015.2020.1867558)

For a bivariate I(1) system, this tests the null of no cointegration (rank `r = 0`).
Alongside the standard Johansen **trace** test, it computes the **adaptive** LR test,
which weights by a nonparametric estimate of the time-varying (spot) error covariance
and so gains power when volatility is nonstationary. Each test is reported with
**wild-bootstrap** and **volatility-bootstrap** p-values. This is a ready-to-use
implementation, not replication code for the paper's tables.

## Requirements

[Ox](https://oxlang.dev) — the free **Ox Console** or **OxMetrics**. No extra packages
(only `oxstd`).

## Usage

```ox
#include <oxstd.oxh>
#include "adaptcoint.ox"

main()
{
    decl X = loadmat("mydata.mat");   // n x 2 matrix: the two I(1) series
    decl r = coint_test(X, 5);        // second argument = max lag for BIC selection
    println("k trace adaptLR p_wb_trace p_wb_adapt p_vb_trace p_vb_adapt");
    println(r);
}
```

`coint_test(X, kmax)` returns a 1×7 row:

| column | meaning |
|--------|---------|
| 1 | `k` — lag order chosen by BIC |
| 2 | trace statistic (Johansen) |
| 3 | adaptive LR statistic |
| 4 | wild-bootstrap p-value, trace |
| 5 | wild-bootstrap p-value, adaptive |
| 6 | volatility-bootstrap p-value, trace |
| 7 | volatility-bootstrap p-value, adaptive |

Reject the no-cointegration null at level α when the p-value < α. The **adaptive test
with the volatility bootstrap** (column 7) is the paper's recommended procedure.

See [`demo.ox`](demo.ox) for a runnable example that simulates a system under the null
and under a cointegrated alternative.

## How to cite

```bibtex
@article{boswijk2022adaptive,
  title   = {Adaptive testing for cointegration with nonstationary volatility},
  author  = {Boswijk, H. Peter and Zu, Yang},
  journal = {Journal of Business and Economic Statistics},
  volume  = {40},
  number  = {2},
  pages   = {744--755},
  year    = {2022},
  doi     = {10.1080/07350015.2020.1867558}
}
```

## Notes

The library (`adaptcoint.ox`) is extracted from the paper's own Monte Carlo code
(`lib2020.ox`, `lib2020k1.ox`, and the driver's `pvals`/`pvalsk1`) and consolidated
into one callable file with a single `coint_test` entry point. It was verified to
reproduce the paper's Monte Carlo output **bit-for-bit** (the size and power values of
Table 1). The current implementation covers the bivariate (`p = 2`) case with restricted
constant, tested for rank `r = 0`. Bootstrap p-values use 499 replications.

## License

MIT — see [LICENSE](LICENSE).
