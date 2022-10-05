options(scipen=999)
library(pscl)
library(ggplot2)
library(tidyverse)
library(metRology)
library(latex2exp)
library(LaplacesDemon)
library(invgamma)
library(TeachingDemos)
my_theme <-theme(
    legend.position='none',
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size = 1)
)
data(absentee)
View(absentee)


# Example 2.14 ====
#prior
n0 = 5
mu0 = 0
n = 21
ybar = -5.8
s = 1161.3
v0 = 6.2
sig2_0 = 47.07

#posterior
mu1 = (n0*mu0 + n*ybar)/(n0+n)
n1 = n0+n
v1 = v0 +n
v1sig2_1 = v0*sig2_0 + s + (n0*n/(n0+n))*(mu0-ybar)^2
sigsq
postvar = (v1sig2_1/2)/(v1/2-1)

E <- function(n0, mu0, n, ybar){
    expectation <- (n0*mu0+n*ybar)/(n0+n)
    expectation
}

n_0 <- seq(0,10000,0.1)
posterior <- E(n_0,0,21,-5.8)  


lower <- double()
upper <- double()
data <- data.frame(n_0, posterior)



ggplot(data,aes(x=n_0, y = posterior)) + 
    geom_line()+ 
    scale_x_continuous(trans='log10', breaks = c(1,10,100,1000,10000))+
    ylim(-8,0)



# Example 2.15 ====

absentee
# start running this
absentee$a <- absentee$absdem/(absentee$absdem+absentee$absrep)
absentee$m <- absentee$machdem/(absentee$machdem+absentee$machrep)
absentee$y <- (2*absentee$a-1)*100
absentee$x <- (2*absentee$m-1)*100
View(absentee)

## Prior densities ====

### prior mean ----
b0 <- c(0,1)
nu0 <- 6.2
sigsq0 <- 47.07

### prior variance

q0 <- 25
qcurl0 <- qt(0.975,df = nu0)
lambda0 <- (q0/qcurl0)^2/sigsq0

q1 <- 1
qcurl1 <- qt(0.99,df = nu0)
lambda1 <- (q1/qcurl1)^2/sigsq0

B0 <- diag(c(lambda0,lambda1),2,2)
# hpd for beta
hpd(qt.scaled, df = nu0, mean = 0, sd = sqrt(sigsq0*lambda0))
hpd(qt.scaled, df = nu0, mean = 1, sd = sqrt(sigsq0*lambda1))
25+1.78*1.45
### prior sig^2 ----

# mean of prior sig^2

expect <- function(x){x*dinvgamma(x,nu0/2, nu0*sigsq0/2)}
integrate(expect,0,Inf)

temp <- c()
for (i in 1:1000){
    temp[i] <- mean(rinvgamma(2*i-1,nu0/2,nu0*sigsq0/2))
}

plot(temp, type = 'line')

#confirming: mean = 69.48429
prior.meansigsq <- (nu0*sigsq0/2)/(nu0/2-1)
prior.modesigsq <- nu0*sigsq0/(nu0+2)

# confidence interval for prior sigsq
prior.sigsq.hpd <- hpd(qinvgamma, shape = nu0/2, rate = nu0*sigsq0/2) # this is consistent with the paper


## PLOT: MLE vs prior ----

mlevsprior <- ggplot(absentee, aes(x,y)) + 
    
    xlab("Democratic margin in machine ballots (x)") +
    ylab("Democratic margin in absentee ballots (y)") +
    
    geom_point(shape = 1, size = 3) + 
        geom_point(data = absentee[22,], size = 3) + 
            geom_text(data = absentee[22,], 
                      aes(label = "Disputed Election"),
                      hjust=-0.1,vjust=0.4, size = 3.5)+
    
    stat_function(fun = function(x){x}, 
                  aes(colour = "Prior Mean")) +
    geom_smooth(data = absentee[1:21,],method = "lm", se = FALSE,
                  aes(colour = "MLE"))+
    scale_colour_manual("",values = c("Prior Mean" = "black", "MLE" = "blue"))
                        


## Data and Likelihood ====

n = 21

X <- matrix(c(rep(1,21),absentee$x[1:21]), ncol = 2)
y <- absentee$y[1:21]

betahat <- solve(t(X)%*%X)%*%t(X)%*%y
solve(t(X)%*%X)

# sum of squares residuals
S <- t(y - X%*%betahat)%*%(y - X%*%betahat) %>%
    as.numeric()


# MLE of sig squared
sigsqMLE <- S/n

# standard error of betahat
varcovar<- sigsqMLE*solve(t(X)%*%X)
sebeta_0 <- sqrt(varcovar[1,1])
sebeta_1 <- sqrt(varcovar[2,2])


## Posterior densities ----

b1 <- solve((solve(B0)+t(X)%*%X))%*%(solve(B0)%*%b0+t(X)%*%X%*%betahat)
B1 <- solve((solve(B0)+t(X)%*%X))

nu1 <- nu0 + n

r <- t(b0-betahat)%*%solve(B0+solve(t(X)%*%X))%*%(b0-betahat) %>%
    as.numeric()

nu1sigsq1 <- nu0*sigsq0 + S + r
sigsq1 <- nu1sigsq1/nu1

# mean of posterior sig^2
post.meansigsq <- (nu1sigsq1/2)/(nu1/2-1)
post.modesigsq <- nu1sigsq1/(nu1+2)

# confidence interval for posterior sigsq
post.sigsq.hpd <- hpd(qinvgamma, shape = nu1/2, rate = nu1sigsq1/2) # this is consistent with the paper



## Predictions for the disputed election ----
xcurl <- c(1,-1.45)
prior.ycurl <- b0%*%xcurl
post.ycurl  <- t(b1)%*%xcurl




## Figure 2.15 graphs ----


#1 marginal prior of mean

priorbeta <- data.frame(beta0 = rt.scaled(1000, df = nu0, mean = b0[1], sd = sqrt(sigsq0*lambda0)),
                beta1 = rt.scaled(1000, df = nu0, mean = b0[2], sd = sqrt(sigsq0*lambda1)))

ggplot(priorbeta, aes(x = beta0, y = beta1)) + my_theme +
    xlim(-25,25) + ylim(0,2)+ 
    labs(title = "Prior", x = TeX("$\\beta_0$"), y = TeX("$\\beta_1$"))+
    
    stat_density_2d(aes(fill = ..density..), 
                    geom = "raster", 
                    contour = FALSE,
                    n=1000) +
    scale_fill_distiller(palette = "Greys", direction = 1) +
    geom_point(data = as.data.frame(t(colMeans(priorbeta))), size = 2, color ='white')
    
# 2 marginal posterior of mean

postbeta<- data.frame(rmvt(n=1000, mu = b1, S= round(sigsq1*B1,4), df = nu1)) %>%
    rename("beta0"= `X1` , "beta1"=`X2` )

ggplot(postbeta, aes(x = beta0, y = beta1)) +  my_theme+
    xlim(-25,25) + ylim(0,2)+ 
    labs(title = "Posterior", x = TeX("$\\beta_0$"), y = TeX("$\\beta_1$"))+
    
    stat_density_2d(aes(fill = ..density..), 
                    geom = "raster", 
                    contour = FALSE,
                    n=1000) +
    scale_fill_distiller(palette = "Greys", direction = 1) +
    geom_point(data = as.data.frame(t(colMeans(postbeta))), size = 1.5, color ='white')

# 3 marginal prior and posterior of sig^2


priorsigsq <- list(nu0/2, nu0*sigsq0/2)
postsigsq  <- list(nu1/2, nu1sigsq1/2)


ggplot(data.frame(x = c(0,400, by = 0.1)), aes(x = x)) + my_theme+
    labs(x = TeX("$\\sigma^2$"))+
    
    stat_function(fun = dinvgamma, args = priorsigsq, size = 1)+
    stat_function(fun = dinvgamma, args = postsigsq, size = 1) + 
    geom_area(stat = "function", fun = dinvgamma, args = postsigsq, fill = "sky blue", xlim = post.sigsq.hpd, alpha = 0.5)+
    
    annotate(geom = "text", label = "Prior", x = 60, y = 0.015)+
    annotate(geom = "text", label = "Posterior", x = 200, y = 0.008)
    
 
# 4 predictions
        # to follow...

## TABLE ----
table <- data.frame("prior"=c(b0,prior.meansigsq,prior.ycurl), "posterior" = c(b1,post.meansigsq,post.ycurl))
rownames(table) <- c("intercept", "slope", "sigma^2", "prediction")
table


## SENSITIVITY ANALYSIS graphs ====

# 1
beta.sim <- bind_rows("prior" = priorbeta,"posterior" = postbeta, .id = 'grp')
model <- lm(y~x, data = absentee[1:21,])
MLE <- data.frame(t(coefficients(model)))
colnames(MLE) <- c("beta0", "beta1")

ggplot(beta.sim, aes(x = beta0, y = beta1)) + my_theme + 
    xlim(-25,25) + ylim(0,2)+ 
    labs(x = TeX("$\\beta_0$"), y = TeX("$\\beta_1$"))+
    
    stat_ellipse(data = beta.sim[beta.sim$'grp'=="prior",],
                 level = 0.95, geom = 'polygon', alpha = 0.2)+
    stat_ellipse(data = beta.sim[beta.sim$'grp'=="posterior",],
                 level = 0.95, geom = 'polygon', alpha = 0.5)+
    
    geom_point(data = aggregate(.~grp, data = beta.sim, mean))+
    geom_text(data = aggregate(.~grp, data = beta.sim, mean),aes(label=paste(grp,"mean")), hjust = -0.1)+
    
    geom_point(data = MLE, color = 'blue',shape = 18)+
    geom_text(data = MLE, aes(label = "MLE"),hjust = 1.2,vjust = 0, color = 'blue')
# 2
# 3/4 
b1.sens <- data.frame()

for (k in seq(0.1,10000,0.1)){
    b1.k <- solve((k*solve(B0)+t(X)%*%X))%*%(k*solve(B0)%*%b0+t(X)%*%X%*%betahat)
    b1.sens <- rbind(b1.sens,c(k,t(b1.k)))
}
colnames(b1.sens) <- c("k","beta0", "beta1")
b1.sens$ycurl <- b1.sens$beta0+absentee$x[22]*b1.sens$beta1

# 3 - beta_0
ggplot(b1.sens, aes(x = k, y = beta0))+ 
    my_theme+ ylim(-5,0) +
    xlab(TeX("$\\kappa"))+ ylab(TeX("$\\beta_0"))+
    geom_line(size=1)+
    scale_x_continuous(trans='log10')

# 4 - beta_1
ggplot(b1.sens, aes(x = k, y = beta1))+ 
    my_theme+ ylim(-5,0) +
    xlab(TeX("$\\kappa"))+ ylab(TeX("$\\beta_1"))+
    geom_line()+
    scale_x_continuous(trans='log10')


# EXERCISE 2.18

## Predictions for the disputed election ----

ggplot(b1.sens, aes(x = k, y = ycurl))+ 
    my_theme+ylim(-6.5,0)+
    xlab(TeX("$\\kappa"))+ ylab(TeX("$\\tilde{y}"))+
    geom_line()+
    scale_x_continuous(trans='log10')+
    
    geom_point(data = b1.sens[b1.sens$k==1,],size = 2)+
    annotate("text", label = "posterior prediction = -6.16",x = 1, y = -5.5)+
    annotate("text", label = TeX("$\\kappa = 1$"),x = 1, y = -5.7)+
    
    geom_hline(yintercept = -1.45, linetype= 'dashed')+
    annotate("text", label = "Prior prediction = -1.45",x = 500, y = -1.6)+
    annotate("text", label = TeX("$\\kappa \\rightarrow \\infty$"),x = 500, y = -1.8) +
    
    geom_hline(yintercept = c(1,-1.45)%*%betahat, linetype = 'dashed')+
    annotate("text", label = TeX("$OLS: \\hat{y} = -6.30$"),x = 0.1, y = -6.08, hjust = 0)+
    annotate("text", label = TeX("$\\kappa \\rightarrow 0$"),x = 0.1, y = -6.4, hjust = 0) +
    
    annotate("text", label = "Disputed election", x = 0.5, y = 0)+
    annotate("text", label = TeX("$y_{22}=58.01$"), x = 0.5, y = -0.25)+
    geom_segment(aes(x = 0.2, y = -0.5, xend = 0.2, yend = 0),
                 arrow = arrow())
    absentee[22,]

c(1,-1.45)%*%betahat
model <- lm(y~x, data = absentee[1:21,])

predict(model, newdata=absentee[22,])
ggplot(absentee, aes(x = x, y = y)) +geom_point()

### biased
biased <-  solve((1000*solve(B0)+t(X)%*%X))%*%(1000*solve(B0)%*%c(40,1)+t(X)%*%X%*%betahat)
biased.B1 <-  solve((1000*solve(B0)+t(X)%*%X))
c(1,-1.45)%*%biased 
xcurl%*%biased.B1
biased.parms <- list(df = nu1, mean = xcurl%*%biased, sd=sqrt(sigsq1*(xcurl%*%biased.B1%*%xcurl+1)))
biased.hpd <- hpd(qt.scaled, df = nu1, mean = xcurl%*%biased, sd=sqrt(sigsq1*(xcurl%*%biased.B1%*%xcurl+1)))

ggplot(data.frame(x = c(0,400, by = 0.1)), aes(x = x)) + my_theme+ xlim(-10,100)+
    labs(x = "Democratic Margin")+
     
    stat_function(fun = dt.scaled, args = biased.parms, size = 1)+
    geom_area(stat = "function", fun = dt.scaled, args = biased.parms, fill = "sky blue", xlim = biased.hpd, alpha = 0.5)+
    geom_vline(xintercept = 58.01)
 


mlevsprior <- ggplot(absentee, aes(x,y)) + 
    
    xlab("Democratic margin in machine ballots (x)") +
    ylab("Democratic margin in absentee ballots (y)") +
    
    geom_point(shape = 1, size = 3) + 
    geom_point(data = absentee[22,], size = 3) + 
    geom_text(data = absentee[22,], 
              aes(label = "Disputed Election"),
              hjust=-0.1,vjust=0.4, size = 3.5)+
    
    stat_function(fun = function(x){x}, 
                  aes(colour = "Prior Mean")) +
    geom_smooth(data = absentee[1:21,],method = "lm", se = FALSE,
                aes(colour = "MLE"))+
    stat_function(fun = function(x){biased[1] + biased[2]*x}, 
                  aes(colour = "New Prior"))+
    scale_colour_manual("",values = c("Prior Mean" = "black", "MLE" = "blue", "New Prior" = 'red'))+
    theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = 1)
    )
