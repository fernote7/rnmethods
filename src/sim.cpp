
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sim(NumericMatrix out, float X0, float mu, float Dt, float sqrdt, float sigma) {
    int nrow = out.nrow(), ncol = out.ncol();



    for (int j = 0; j < ncol; j++) {
        out(0,j) = X0;
        for (int i = 0; i < nrow; i++) {
            out(i+1,j) = out(i,j) + mu * out(i,j) * Dt + sigma * out(i,j) * sqrdt * out(i+1,j);
        }
    }

    return out;
}
