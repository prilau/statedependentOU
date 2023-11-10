a = 0.1234
b = 0.2345
c = 0.3456
theta0 = 0.5
theta1 = 2
halflife = 0.35
alpha = log( 2 ) / halflife
stV = 0.0625
sigma2 = 2 * alpha * stV

#branch A
vL = sigma2 / ( 2 * alpha ) * ( exp(2 * alpha * 0.7) - 1 )
dL = 0
varL = vL + dL * exp(2 * alpha * 0.7)
dL = varL
muL = exp(0.7 * alpha) * (a  - theta0)  + theta0


#branch B 
vR = sigma2 / ( 2 * alpha ) * ( exp(2 * alpha * 0.7) - 1 )
dR = 0
varR = vR + dR * exp(2 * alpha * 0.7)
dR = varR
muR = exp(0.7 * alpha) * (b  - theta0)  + theta0

mu_nAB = ( muL * varR + muR * varL ) / ( varL + varR )
mu_contrast = muL - muR

zL = 1
zR = 1
z_nAB = exp( alpha * 0.7 + alpha * 0.7 ) / ( zL * zR ) * exp( -mu_contrast^2 / ( 2 * ( varL + varR ) ) ) / sqrt( 2 * pi ) / sqrt( varL + varR )
lnl_nAB = log(z_nAB)

#branch AB
vL = sigma2 / ( 2 * alpha ) * ( exp(2 * alpha * 0.3) - 1 )
varL = vL + dL * exp(2 * alpha * 0.3)
dL = varL
muL = exp(0.3 * alpha) * (mu_nAB  - theta0)  + theta0

#branch C
dR = 0
vR = sigma2 / ( 2 * alpha ) * ( exp(2 * alpha * 0.5) - 1 )
varR = vR + dR * exp(2 * alpha * 0.5)
dR = varR
muR = exp(0.5 * alpha) * (c  - theta1)  + theta1
vR = sigma2 / ( 2 * alpha ) * ( exp(2 * alpha * 0.5) - 1 )
varR = vR + dR * exp(2 * alpha * 0.5)
dR = varR
muR = exp(0.5 * alpha) * (muR  - theta0)  + theta0


mu_nABC = ( muL * varR + muR * varL ) / ( varL + varR )
mu_contrast = muL - muR
zL = 1
zR = 1
z_nABC = exp( alpha * 1.0 + alpha * 1.0 ) / ( zL * zR ) * exp( -mu_contrast^2 / ( 2 * ( varL + varR ) ) ) / sqrt( 2 * pi ) / sqrt( varL + varR )
lnl_nABC = log(z_nABC) + lnl_nAB
