/* 
Eugene Cho
Math 2BL
Proffessor: Jeff Anderson
18 November 2025
*/
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

struct MCResult { double price; double stderr; double conf_low; double conf_high; };

MCResult price_european_call_gbm(double S0, double K, double r, double mu,
                                 double sigma, double T, int nPaths,
                                 uint64_t seed = std::random_device{}())
{
    mt19937_64 gen(seed);
    normal_distribution<double> N01(0.0, 1.0);
    const double drift = (mu - 0.5 * sigma * sigma) * T;
    const double vol = sigma * sqrt(T);
    const double disc = exp(-r * T);
    double mean = 0.0, M2 = 0.0;
    for (int i = 1; i <= nPaths; ++i) {
        double Z = N01(gen);
        double ST = S0 * exp(drift + vol * Z);
        double payoff = max(ST - K, 0.0);
        double delta = payoff - mean; mean += delta / i; M2 += delta * (payoff - mean);
    }
    double variance = (nPaths > 1) ? (M2 / (nPaths - 1)) : 0.0;
    double stddev = sqrt(variance);
    double stderr = stddev / sqrt((double)nPaths);
    double price = disc * mean;
    double ci = 1.96 * disc * stderr;
    return { price, disc * stderr, price - ci, price + ci };
}

// Price a 50/50 basket call on two assets with correlation rho (uses 2x2 Cholesky-like factor)
MCResult price_basket_call_correlated(double S0a, double S0b, double K, double r,
                                     double muA, double muB, double sigmaA, double sigmaB,
                                     double rho, double T, int nPaths,
                                     uint64_t seed = std::random_device{}())
{
    mt19937_64 gen(seed); normal_distribution<double> N01(0.0,1.0);
    const double driftA = (muA - 0.5*sigmaA*sigmaA)*T, driftB = (muB - 0.5*sigmaB*sigmaB)*T;
    const double volA = sigmaA*sqrt(T), volB = sigmaB*sqrt(T), disc = exp(-r*T);
    const double L11 = 1.0, L21 = rho, L22 = sqrt(max(0.0, 1.0 - rho*rho)); // L * z produces correlated normals
    double mean = 0.0, M2 = 0.0;
    for (int i = 1; i <= nPaths; ++i) {
        double z1 = N01(gen), z2 = N01(gen);
        double c1 = L11 * z1;
        double c2 = L21 * z1 + L22 * z2;                 // 2x2 multiply (L · z)
        double STa = S0a * exp(driftA + volA * c1);
        double STb = S0b * exp(driftB + volB * c2);
        double basket = 0.5 * STa + 0.5 * STb;
        double payoff = max(basket - K, 0.0);
        double delta = payoff - mean; mean += delta / i; M2 += delta * (payoff - mean);
    }
    double variance = (nPaths > 1) ? (M2 / (nPaths - 1)) : 0.0;
    double stddev = sqrt(variance), stderr = stddev / sqrt((double)nPaths);
    double price = disc * mean, ci = 1.96 * disc * stderr;
    return { price, disc * stderr, price - ci, price + ci };
}

int main() {
    double S0 = 100, K = 100, r = 0.05, mu = 0.05, sigma = 0.2, T = 1.0;
    int nPaths = 200000;
    auto res = price_european_call_gbm(S0, K, r, mu, sigma, T, nPaths);
    auto res2 = price_basket_call_correlated(100,100,100,r,mu,mu,0.2,0.2,0.6,T,nPaths);
    cout.setf(ios::fixed); cout<<setprecision(6);
    cout << "MC price (single asset): " << res.price << "\n"
         << "MC price (50/50 correlated basket): " << res2.price << "\n";
         << "Standard Error  : " << res.stderr << "\n"
         << "95% Confidence Interval  : [" << res.conf_low << ", " << res.conf_high << "]\n";
    return 0;
}
