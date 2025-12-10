# Monte Carlo Simulation with Geometric Brownian Motion

The modeling project that I worked on was a Monte Carlo simulation Geometric Brownian Motion. I used C++ to code this model. 
The model prices a European call option at a fixed price basically by asking “if the future stock price could move in different random ways,
what would the average payoff be?” and gives that price as the Monte Carlo price. By simulating a number of future outcomes, the program 
finds the differences between the stock price and the strike price and averages the differences. Because the model calculates using options 
which are only rights to purchase rather than actual stocks themselves, if the stock ends up below the strike price, the difference, or 
Monte Carlo price, of that instance is set to 0. 

I used vector-matrix multiplication between a 2x2 lower-triangular matrix, L and a vector z, consisting of independent normals. The vector
-matrix multiplication turns the independent normal sample into a correlated sample which allows me to determine how the assets move together.

## What is output: 
  MC price (single asset)<br>
  MC price (50/50 correlated basket)<br>
  Standard Error<br>
  95% Confidence Interval<br>

## Math Model
  Geometric Brownian Motion:<br>
  S_T = S_0 * exp( (mu - 0.5*sigma^2)*T + sigma*sqrt(T)*Z )<br>
  	•	S0 = initial price<br>
	•	mu = drift<br>
	•	sigma = volatility<br>
	•	T = time to maturity<br>
	•	Z ~ N(0, 1)<br>
  Asset prices face daily shocks which are normally distributed with large jumps being quite rare. We generate the normal distribution of 
  shocks:<br>
  	normal_distribution<double> N01(0.0, 1.0);<br>
  A mean of 0 with a standard deviation of 1. Geometric Brownian Motion assumes this normal distribution of change for an asset over time.
  By simulating the shocks nPaths amount of times for nPaths amount of futures, we can avearge the Monte Carlo prices along the Geometric
  Brownian Motion.

## Initial Values

Hypothetical initial asset and strike values for clearest display of volatility:<br>
  S0 = 100			Initial Price<br>
  K = 100			Strike Value<br>
  
  Risk neutral setting:<br>
  r = 0.05			Risk Free Interest Rate<br>
  mu = 0.05			Drift<br>

  sigma = 0.2    	Volatility<br>
  T = 1.0			Time to Maturity<br>
  nPaths = 200000	Number of Paths<br>
  
### European Call Price:
  max(S_T - K, 0)<br>

### Discounted Price
  exp(-r*T) * E[payoff]<br>

## Correlated Basket Option

Because of doubled diversification, the correlated 50/50 basket options have less volatility. This means that with lower volatility,
there are lower chances of ending above the strike price, thus lower prices. We should see this expressed in the output of the model.

### Two correlated standard normals are generated using a 2x2 factorization:<br>
  c1 = z1<br>
  c2 = rho*z1 + sqrt(1 - rho^2)*z2<br>

### Asset Paths
  S_T^A = S0A * exp(driftA + volA * c1)<br>
  S_T^B = S0B * exp(driftB + volB * c2)<br>

### Basket Payoff
  max( 0.5*S_T^A + 0.5*S_T^B - K, 0 )<br>

### Compile
  g++ -std=c++17 -O2 monteCarloGBM.cpp -o monte_carlo<br>
### Run
  ./monte_carlo<br>

### Expected output:
  MC price (single asset): 10.442752<br>
  MC price (50/50 correlated basket): 9.673678<br>
  Standard Error  : 0.032894<br>
  95% Confidence Interval  : [10.378280, 10.507224]<br>

## Verification

### Single asset verification
  The MC price of the single asset can be verified through calculating and comparing the same MC price with the Black-Scholes formula.
  With same variables held constant:<br>
  double S0 = 100, K = 100, r = 0.05, mu = 0.05, sigma = 0.2, T = 1.0;<br>
  int nPaths = 200000;<br>
  The Black-Scholes formula will yield C ~= 10.4492<br>
  Which our expected output will contain in its confidence interval.<br>

The following is a handwritten verification of this value. 
[blackscholesverification.pdf](https://github.com/user-attachments/files/24065755/blackscholesverification.pdf)



    
### Convergence verification
  Monte carlo simulations should stabilize as the number of paths increases. i.e. Law of large numbers. To verify this, We can change
  the values of nPaths in increasing increments to see how the prices tend towards 10.45

#### Some resources I used:

https://www.notion.so/Sources-29baad7fb979800caf64dcab92dac7c9?source=copy_link
  
