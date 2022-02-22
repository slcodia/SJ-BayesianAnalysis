# EXERCISE 1.2
# Contour plot for showing posterior probabilities

library(plotly)

# PRELIMINARY
prior <- 0.03

FP <- c(seq(0.01, 0.30, by=0.01))
FN <- c(seq(0.01, 0.30, by=0.01))



# FIRST TEST =========================================
post1 <- matrix(data = NA,length(FP),length(FP))

# computation
for (i in 1:length(FP)){
    for (j in 1:length(FP)){
        post1[i,j] <- prior*(1-FN[i])/(prior*(1-FN[i])+(1-prior)*FP[j])
    }
}
# plot
fig1 <- plot_ly(
              x = FP,
              y = FN,
              z = post1, 
              type = "contour" 
)%>%
    layout(title = "Probability of substance use after drug test",
           xaxis = list(title = "False Positive Rate"),
           yaxis = list(title = "False Negative Rate")
)%>%
    colorbar(title = list(text = "posterior \nprobability",  font = list(size = 10)), len = 1)

fig1

# SECOND TEST ================================================================================

post2 <- matrix(,length(FP),length(FP))

# computation: posterior from first test becomes prior in the second test
for (i in 1:length(FP)){
    for (j in 1:length(FP)){
        post2[i,j] <- post1[i,j]*(1-FN[i])/(post1[i,j]*(1-FN[i])+(1-post1[i,j])*FP[j])
    }
}

fig2 <- plot_ly(
    x = FP,
    y = FN,
    z = post2, 
    type = "contour" 
)%>%
    layout(title = "Probability of substance use after SECOND drug test",
           xaxis = list(title = "False Positive Rate"),
           yaxis = list(title = "False Negative Rate")
    )%>%
    colorbar(title = list(text = "posterior \nprobability",  font = list(size = 10)), len = 1)

fig2

# THIRD TEST ================================================================================

post3 <- matrix(,length(FP),length(FP))

# computation: posterior from second test becomes prior in the third test
for (i in 1:length(FP)){
    for (j in 1:length(FP)){
        post3[i,j] <- post2[i,j]*(1-FN[i])/(post2[i,j]*(1-FN[i])+(1-post2[i,j])*FP[j])
    }
}

fig3 <- plot_ly(
    x = FP,
    y = FN,
    z = post3, 
    type = "contour" 
)%>%
    layout(title = "Probability of substance use after THIRD drug test",
           xaxis = list(title = "False Positive Rate"),
           yaxis = list(title = "False Negative Rate")
    )%>%
    colorbar(title = list(text = "posterior \nprobability",  font = list(size = 10)), len = 1)

fig3

