############################################################
### First EBola model, trying to see he effects of sGP #####
############################################################

###################################
# FUNCTION DEFINITIONS
###################################

###
# ebola(t,x,parms)
# Use:	Function to calculate derivative of equations in the basic HIV model.
#       Note that this is a variant with no explicit variable for the virus level.
# Input:
# 	t: time (not used here, because there is no explicit time dependence)
# 	x: vector of the current values of all variables
# 	parms: 	dummy variable, which is not used here (normally used to pass on
#		parameter values, but not needed here because parameters are defined globally)
# Output:
#	der: vector of derivatives

ebola <- function(t,x,parms){
  with(as.list(c(x)),{# this structure allows us to refer to the elements of x directly
    dT <- sigma - d*T - b*T*V
    dI <- b*T*V - d_I*I -k*E*I
    dV <- p*I -d_V*V
    dE <- alpha*E*I -d_E*E
    der <- c(dT,dI,dV,dE)
    list(der)
  })
}

n.integrate <- function(start,end,intstep=10,init.x,model){
  t.out <- seq(start,end,length=(end-start)*intstep+1)
  as.data.frame(lsoda(init.x,t.out,model,parms=parms))
}

###########################
# MAIN PROGRAM
###########################

### LOAD LIBRARIES
#load R library for ordinary differential equation solvers
if (!require("deSolve"))
{
  install.packages("deSolve", dependencies = TRUE)
  library(deSolve)
}

### INITIALIZE PARAMETER SETTINGS
## Model parameters (defined on the time scale of days)
#  I define parameters globally to be able to change them between runs easily

sigma=10 #all params are per day
d=0.01
d_I = 0.1
d_E = 0.1
d_V = 3
p = 100
b=0.01
k = 1
alpha = 2

parms <- c() # dummy-variable needed for lsoda

# Initial conditions (i.e. the values of all variables at time = start)
I0 <- 0   # to model initial infection
T0 <- sigma/d   # to start from uninfected steady state
E0 <- 0.1
V0 <- 1
init.x <- c(T=T0, I=I0, V=V0, E=E0)

####################################################
### PART I
### Simulation of initial infection

simlength <- 100
out <- n.integrate(start=0,end=simlength,init.x=init.x,model=ebola)
attach(out)  # this allows us to refer to data columns in out directly

# Plot time course
plot(c(0,simlength),c(0,T0),type="n",xlab="time (days)",ylab="", ylim = c(0,2000))
lines(time,rep(0,length(time)),col="grey",lwd=1,lty=2) # This command draws a line with f(t)=0 to better observe if a population converges against 0, i.e. if the population will go extinct
lines(time,T,lwd=2)
lines(time,I,col="red",lwd=2)
lines(time,V,col="green",lwd=2)
lines(time,E,col="blue",lwd=2)
legend("topright",c("T","I", "V", "E"),lwd=c(2,2,2,2),col=c("black","red", "green", "blue"))
