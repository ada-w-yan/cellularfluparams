# the single-strain TEIV model with stages and dual measurement

# number of stages
n_L <- user()
n_I <- user()

# initial conditions
initial(T1) <- T_0
initial(L[1:n_L]) <- L_0[i]
initial(I[1:n_I]) <- I_0[i]
initial(V_inf) <- V_inf_0
initial(V_tot) <- V_tot_0

# equations
deriv(T1) <- - max(0, beta * T1 * V_inf)
deriv(L[1]) <- max(0, beta * T1 * V_inf) - max(0, k1 * n_L * L[i])
deriv(L[2:n_L]) <- max(0, k1 * n_L * L[i-1]) - max(0, k1 * n_L * L[i])
deriv(I[1]) <- max(0, k1 * n_L * L[n_L]) - max(0, delta * n_I * I[i])
deriv(I[2:n_I]) <- max(0, delta * n_I * I[i-1]) - max(0, delta * n_I * I[i])
deriv(V_inf) <- max(0, p_inf * sum(I)) - max(0, c_inf * V_inf) - max(0, beta_inf * T1 * V_inf)
deriv(V_tot) <- max(0, p_tot * sum(I)) - max(0, c_tot * V_tot) - max(0, beta_tot * T1 * V_inf)

dim(L) <- n_L
dim(I) <- n_I
dim(L_0) <- n_L
dim(I_0) <- n_I

# parameter values
beta <- user()
beta_inf <- user()
beta_tot <- user()
k1 <- user()
delta <- user()
p_inf <- user()
p_tot <- user()
c_inf <- user()
c_tot <- user()
T_0 <- user()
V_inf_0 <- user()
V_tot_0 <- user()
L_0[] <- user()
I_0[] <- user()
