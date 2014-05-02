# New to old ratio

### need to get the neg and pos working



# download pca object and count table
# make sure that after size selection we get the right amount

# L2 + MIR 
# anything species specific 
# in human this will probably be alu

setwd("~/Desktop/new_old_ratio/")

mb <- 1000000
bin.size = 500000

dat <- "Human/AllBinCounts.txt"
count <- read.table(dat, header = TRUE)
count <- count[count$Known >= bin.size,]
countmb <-  count
countmb[,5:(length(countmb)-1)] <- countmb[,5:(length(countmb)-1)]/mb   

load("Human/pca")

dim(cpca$x)[1] == dim(count)[1]


plot(cpca$x[,2], ((count$SINE2_MIR/mb) + (count$LINE_L2/mb)) / (count$SINE1_7SL/mb))

plot(cpca$x[,2], ((countmb$SINE2_MIR + countmb$LINE_L2)))
cor(cpca$x[,2], ((countmb$SINE2_MIR + countmb$LINE_L2)))

plot(cpca$x[,2], countmb$SINE1_7SL)



biplot(cpca)





# how to make the wigle plot
# calculate the 95% conf interval on chr 3 for the pc2 score with a running ave of 5 bins
# then repeat for cow



pc2.chr3 <- data.frame(countmb, cpca$x[,2]) 
pc2.chr3 <- pc2.chr3[pc2.chr3$Chromosome == "chr1",]
pc2 <- pc2.chr3[,length(pc2.chr3)]
# because of the way it turned out we multiply by negative 1
if(cor(cpca$x[,2], ((countmb$SINE2_MIR + countmb$LINE_L2))) < 0){pc2 <- pc2 *-1}

# get a bin and get the bins next to it and calculate the mean score and the sd and thereby get the confidence interval
# iitialize a dataframe to add things to
# 


zero <- rep(0, length(pc2))
wiggle <- data.frame(pc2_score = zero, mean = zero , bottom = zero , top = zero)

for(i in seq(along=pc2)){
	if(i == 1){
		a <- mean(pc2[i:(i+1)])
		s <- sd(pc2[i:(i+1)])
		n <- 2
		error <- qnorm(0.975)*s/sqrt(n)
		left <- a-error
		right <- a+error
	}else if(i == 2){
		a <- mean(pc2[(i-1):(i+1)])
		s <- sd(pc2[(i-1):(i+1)])
		n <- 3
		error <- qnorm(0.975)*s/sqrt(n)
		left <- a-error
		right <- a+error
	}else if(i == length(pc2) - 1){
		a <- mean(pc2[(i-1):(i+1)])
		s <- sd(pc2[(i-1):(i+1)])
		n <- 3
		error <- qnorm(0.975)*s/sqrt(n)
		left <- a-error
		right <- a+error
	}else if(i == length(pc2)){
		a <- mean(pc2[(i-1):i])
		s <- sd(pc2[(i-1):i])
		n <- 2
		error <- qnorm(0.975)*s/sqrt(n)
		left <- a-error
		right <- a+error
	}else{
		a <- mean(pc2[(i-2):(i+2)])
		s <- sd(pc2[(i-2):(i+2)])
		n <- 5
		error <- qnorm(0.975)*s/sqrt(n)
		left <- a-error
		right <- a+error
	}
	
	wiggle$pc2_score[i] <- pc2[i]
	wiggle$mean[i] <- a 
	wiggle$bottom[i] <- left		
	wiggle$top[i] <- right
		
 	
}


plot(wiggle$mean, type = "l")
lines(wiggle$top, col = 2)
lines(wiggle$bottom, col = 2)
lines(zero, col = 3)
lines(wiggle$pc2_score)

lines(scale((pc2.chr3$SINE2_MIR) + (pc2.chr3$LINE_L2))[,1], col = 5, lwd = 3)
lines(scale(pc2.chr3$SINE1_7SL)[,1], col = 6, lwd = 3)

lines(scale(pc2.chr3$SINE1_7SL / (pc2.chr3$SINE2_MIR + pc2.chr3$LINE_L2))[,1], col = 4, lwd = 3)
lines(((pc2.chr3$SINE2_MIR + pc2.chr3$LINE_L2) / pc2.chr3$SINE1_7SL), col = 4, lwd = 3)


lines(scale((pc2.chr3$SINE2_MIR + pc2.chr3$LINE_L2) - pc2.chr3$SINE1_7SL)[,1], col = 4, lwd = 3)

