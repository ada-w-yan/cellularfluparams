# the mock yield model

# initial conditions
initial(V) <- V_0

# equations
deriv(V) <- - c * V

# parameter values
V_0 <- user()
c <- user()
