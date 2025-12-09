# Monte Carlo Simulation with Geometric Brownian Motion

The modeling project that I worked on was a Monte Carlo simulation Geometric Brownian Motion. I used C++ to code this model. 
The model prices a European call option at a fixed price basically by asking “if the future stock price could move in different random ways,
what would the average payoff be?” and gives that price as the Monte Carlo price. By simulating a number of future outcomes, the program 
finds the differences between the stock price and the strike price and averages the differences. Because the model calculates using options 
which are only rights to purchase rather than actual stocks themselves, if the stock ends up below the strike price, the difference, or 
Monte Carlo price, of that instance is set to 0. 

I used vector-matrix multiplication between a 2x2 lower-triangular matrix, L and a vector z, consisting of independent normals. The vector-matrix multiplication turns the independent normal sample into a correlated sample which allows me to determine how the 
assets move together.

## What is output: 
  MC price (single asset)
  MC price (50/50 correlated basket)
  Standard Error
  95% Confidence Interval

## Math Model
  Geometric Brownian Motion:
  S_T = S_0 * exp( (mu - 0.5*sigma^2)*T + sigma*sqrt(T)*Z )
  •	S0 = initial price
	•	mu = drift
	•	sigma = volatility
	•	T = time to maturity
	•	Z ~ N(0, 1)
  
### European Call Price:
  max(S_T - K, 0)

### Discounted Price
  exp(-r*T) * E[payoff]

## Correlated Basket Option

### Two correlated standard normals are generated using a 2x2 factorization:
  c1 = z1
  c2 = rho*z1 + sqrt(1 - rho^2)*z2

### Asset Paths
  S_T^A = S0A * exp(driftA + volA * c1)
  S_T^B = S0B * exp(driftB + volB * c2)

### Basket Payoff
  max( 0.5*S_T^A + 0.5*S_T^B - K, 0 )

### Compile
  g++ -std=c++17 -O2 monte_carlo.cpp -o monte_carlo
### Run
  ./monte_carlo

### Expected output:
  MC price (single asset): 10.442752
  MC price (50/50 correlated basket): 9.673678
  Standard Error  : 0.032894
  95% Confidence Interval  : [10.378280, 10.507224]

## Verification

### Single asset verification
  The MC price of the single asset can be verified through calculating and comparing the same MC price with the Black-Scholes formula.
  With same variables held constant:
  double S0 = 100, K = 100, r = 0.05, mu = 0.05, sigma = 0.2, T = 1.0;
  int nPaths = 200000;
  The Black-Scholes formula will yield C ~= 10.4506
  Which is within 1 standard deviation with our expected output.
  
### Convergence verification
  Monte carlo simulations should stabilize as the number of paths increases. i.e. Law of large numbers. To verify this, We can change
  the values of nPaths in increasing increments to see how the prices tend towards 10.45
  
  
