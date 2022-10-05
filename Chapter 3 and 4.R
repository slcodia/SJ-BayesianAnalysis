# Probset 3
library(ggplot2)
library(patchwork)
library(tidyverse)
options(scipen=999)
install.packages("patchwork")
# 3.3

data <- matrix(c(1,2,7,74),2,2)
rownames(data) <- c("defeated", "not defeated")
colnames(data) <- c("revolution", "no revolution")
data

# revolution | not defeated
theta0 <- data[2,1]/sum(data[2,])
# revolution | defeated
theta1 <- data[1,1]/sum(data[1,])
theta1/theta0

n <- 1000
theta_0 <- rbeta(n,3,75)
theta_1 <- rbeta(n,2,8)

#1 ----
ratio <- theta_1/theta_0


#2 ----
OR <- (theta_1/(1-theta_1))/(theta_0/(1-theta_0))


#3 ----
logOR <- log(OR)

# PLOTS


ratio.plot <- ggplot(data.frame(ratio), aes(ratio))  + ylab("count")+xlim(0,15) +
    geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(ratio)), linetype = "dashed")+
    geom_vline(aes(xintercept = 1))+
    annotate(geom = "text", label = paste0("mean = ", round(mean(ratio),2)), 
             x = mean(ratio)+1, y = 10, hjust = 0)

OR.plot <- ggplot(data.frame(OR), aes(OR))  + ylab("count")+ xlim(0,15) +
    geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(OR)), linetype = "dashed")+
    geom_vline(aes(xintercept = 1))+
    annotate(geom = "text", label = paste0("mean = ", round(mean(OR),2)), 
             x = mean(OR)-3, y = 10)

logOR.plot <- ggplot(data.frame(logOR), aes(logOR)) + ylab("count") + xlim(-5,5)+
    geom_histogram(binwidth=0.1) + 
    geom_vline(aes(xintercept=mean(logOR)), linetype = "dashed")+
    geom_vline(aes(xintercept = 0))+
    annotate(geom = "text", label = paste0("mean = ", round(mean(logOR),2)), 
             x = -2.5, y = 15)

ratio.plot + OR.plot + logOR.plot

    ######


data2 <- data.frame(ratio,OR,logOR) %>% 
    pivot_longer(cols = c(ratio, OR, logOR))
View(data2)

ggplot(data2, aes(x=value)) + xlim(-5,5) +
    geom_histogram(data = subset(data2,name == 'ratio'), fill = "red", alpha = 0.2, binwidth = 0.1) + 
    geom_histogram(data = subset(data2,name == 'OR'), fill = "blue", alpha = 0.2, binwidth = 0.1) +
    geom_histogram(data = subset(data2,name == 'logOR'), fill = "green", alpha = 0.2,  binwidth = 0.1)

data3 <-data.frame(theta_0, theta_1) %>%
    pivot_longer(cols = c(theta_0, theta_1))
ggplot(data3, aes(x=value)) + xlim(0,1) +
    geom_histogram(data = subset(data3,name == 'theta_0'), fill = "red", alpha = 0.2, binwidth = 0.01) + 
    geom_histogram(data = subset(data3,name == 'theta_1'), fill = "blue", alpha = 0.2, binwidth = 0.01)
data3
hist(theta_0)
hist(theta_1)



# 4.12 ################


M = 5000; rho = 0.95; t = 10^5
zbar <- c();ybar <- c()
for (m in 1:M){
    # a: AR(1) process
    z_t <- arima.sim(n = t, list(ar = rho), sd = sqrt(1-rho^2))
    zbar[m] <- mean(z_t)
    # b: independence sampler
    y_t <- rnorm(t)
    ybar[m] <- mean(y_t)

}

monte.data <- data.frame(index = seq(1:M), zbar,ybar) %>%
    pivot_longer(cols = -index, names_to = "series")
  
ggplot(monte.data, aes(x = index, y = value, color = series)) + ylim(-0.06, 0.06) +
    geom_point(alpha = 0.5)

?pivot_longer
var(zbar)/var(ybar)
(1+rho)/(1-rho)


#computing pi

u <- runif(100)
v <- runif(100)

