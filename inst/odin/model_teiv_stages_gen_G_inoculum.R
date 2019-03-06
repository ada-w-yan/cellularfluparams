# the single-strain TEIV model with many latent and infectious stages
# counting generations: g = G > 2

n_L <- user()
n_I <- user()

# equations
deriv(T1) <- -beta * T1 * sum(V)
deriv(L[1, 1:G]) <- beta * T1 * V[j] - k1 * n_L * L[1, j]
deriv(L[2:n_L, 1:G]) <- k1 * n_L * (L[i-1, j] - L[i, j])
deriv(I[1, 1:G]) <- k1  * n_L * L[n_L, j] - delta * n_I * I[1, j]
deriv(I[2:n_I, 1:G]) <- delta * n_I * (I[i-1, j] - I[i, j])
deriv(V[1]) <- - (c_inf + beta_inf * T1) * V[1]
deriv(V[2:(G-1)]) <- p_inf * sum(I[,i-1]) - (c_inf + beta_inf * T1) * V[i]
deriv(V[G]) <- p_inf * (sum(I[,G-1]) + sum(I[,G])) - (c_inf + beta_inf * T1) * V[G]
deriv(W[1]) <- 0
deriv(W[2:(G-1)]) <- p_inf * sum(I[,i-1])
deriv(W[G]) <- p_inf * (sum(I[,G-1]) + sum(I[,G]))

dim(L) <- c(n_L, G)
dim(I) <- c(n_I, G)
dim(V) <- G
dim(W) <- G
dim(L_0) <- c(n_L, G)
dim(I_0) <- c(n_I, G)

# initial conditions
initial(T1) <- T_0
initial(L[1:n_L,1:G]) <- L_0[i,j]
initial(I[1:n_I,1:G]) <- I_0[i,j]
initial(V[1]) <- V_inf_0
initial(V[2:G]) <- 0
initial(W[1]) <- V_inf_0
initial(W[2:G]) <- 0

# parameter values
beta <- user()
beta_inf <- user()
delta <- user()
p_inf <- user()
c_inf <- user()
T_0 <- user()
V_inf_0 <- user()
G <- user()
k1 <- user()
L_0[,] <- user()
I_0[,] <- user()