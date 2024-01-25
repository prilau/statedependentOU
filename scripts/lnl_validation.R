library(ape)
tree <- read.tree("data/1_validation/dummy_threetaxon_stateless.tre")
plot(tree)

a = 0.1234
b = 0.2345
c = 0.3456
alpha = 0.1
theta = 0.1
sigma2 = 0.1

Ya = theta * ( 1 - exp(alpha * 0.7)) + exp( alpha * 0.7 ) * a
Yb = theta * ( 1 - exp(alpha * 0.7)) + exp( alpha * 0.7 ) * b
Va = sigma2 / 2 / alpha * ( exp( 2 * alpha * 0.7) - 1 ) + exp( alpha * 0.7 ) * 0
Vb = sigma2 / 2 / alpha * ( exp( 2 * alpha * 0.7) - 1 ) + exp( alpha * 0.7 ) * 0
contrast_ab = Ya - Yb
VarC_ab = Va + Vb
lnl_ab = log(exp( alpha * 0.7 + alpha * 0.7 - contrast_ab^2 / 2 / VarC_ab ) / sqrt( 2 * pi * VarC_ab ))

ab = ( Ya * Vb + Yb * Va ) / ( Va + Vb )
V_pre_ab = ( Va * Vb ) / ( Va + Vb )

Yc = theta * ( 1 - exp( alpha * 1.0 )) + exp( alpha * 1.0 ) * c
Yab = theta * ( 1 - exp( alpha * 0.3 )) + exp( alpha * 0.3 ) * ab
Vc = sigma2 / 2 / alpha * ( exp( 2 * alpha * 1.0 ) - 1 ) + exp( alpha * 1.0 ) * 0
Vab = sigma2 / 2 / alpha * ( exp( 2 * alpha * 0.3) - 1 ) + exp( alpha * 0.3 ) * V_pre_ab
contrast_abc = Yc - Yab
VarC_abc = Vc + Vab
lnl_ab_c = log(exp( alpha * 1.0 + alpha * 0.3 - contrast_abc^2 / 2 / VarC_abc ) / sqrt( 2 * pi * VarC_abc ))

lnl_a_b_c = lnl_ab + lnl_ab_c

root = ( Yc * Vab + Yc * Vab ) / ( Vc + Vab )
Vroot = ( Vc * Vab ) / ( Vc + Vab )

lnproot = log(exp(( root - theta )^2 / ( 2 * Vroot )) / sqrt( 2 * pi * Vroot ))
lnl_everything = proot + lnl_a_b_c
