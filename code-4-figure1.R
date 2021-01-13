
X1=scan()
# 1: -4  0
# 3: -3  0
# 5: -2  0
# 7: -1  0
# 9: 0  0
# 11: 1  0
# 13: 2 -5
# 15: 3  5
# 17: 12  1
# 19: 
#   Read 18 items
 X_1=matrix(X1,ncol=2,byrow=T)

 par(mfrow=c(1,1))
 
# plot(X_1[,1], X_1[,2], pch=19)
# abline(-1.1276501, 0.3010287, lty=1 ) #deepest RD line from rdepthmedian

 plot(X_1, type="p", pch=16, cex=1, xlab="x-axis", ylab='y-axis', 
            xlim=c(-5,15), ylim=c(-5,15))
 title(main="Regression with no point contaminated")
 lines(X_1[,1],lm(X_1[,2]~X_1[,1])$fitted.value, pch=18, col="red", lty=1, lwd=1)
 abline(h=0,lty=2,pch=19, col="blue")
 abline(h=0,lty=3,pch=20, col="black")
 legend(-5, 15, legend=c('LS',expression("T*"["RD"]),expression("T*"["PRD"])),col=c("red", "blue", "black"),
                lty=1:3, cex=0.7, title="Line types", text.font=4, bg='lightblue')
 text(11,3, "LS", col="red", cex=.8)
 text(-1, -1.5, expression("T*"["RD"]), col="blue", cex=0.8)
 text(8,-1.5, expression("T*"["PRD"]), col="black", cex=0.8)

# > X2=scan()
# 1: -4  0
# 3: -3  0
# 5: -2  0
# 7: -1  0
# 9: 1  0
# 11: 2 -5
# 13: 3  5
# 15: 12  12
# 17: 
#   Read 16 item
X_2=matrix(X2, ncol=2, byrow=T)
plot(X_2, type="p", pch=16, cex=1, xlab="x-axis", ylab='y-axis', 
            xlim=c(-5,15), ylim=c(-5,15))
title(main="Regression with one point contaminated")
lines(X_2[,1],lm(X_2[,2]~X_2[,1])$fitted.value, pch=18, col="red", lty=1, lwd=1)
 abline(h=0,lty=2,pch=19, col="blue")
 abline(h=0,lty=3,pch=20, col="black")
 legend(-5, 15, legend=c("LS",expression("T*"["RD"]),expression("T*"["PRD"])),col=c("red", "blue", "black"),
             lty=1:3, cex=0.7, title="Line types", text.font=4, bg='lightblue')
 text(11,6, "LS", col="red", cex=.8)
 text(-1, -1.5, expression("T*"["RD"]), col="blue", cex=0.8)
 text(8,-1.5, expression("T*"["PRD"]), col="black", cex=0.8)
 
 X3=scan()
 # 1: -4  0
 # 3: -3  0
 # 5: -2  0
 # 7: 3  4.5
 # 9: 3  4
 # 11: 3  3.5
 # 13: 2 -5
 # 15: 3  5
 # 17: 12  1
 # 19: 
 #   Read 18 items
 X_3=matrix(X3, ncol=2, byrow=T)
 plot(X_3, type="p", pch=16, cex=1, xlab="x-axis", ylab='y-axis', xlim=c(-5,15), 
      ylim=c(-5,15))
 title(main="Regression with three points contaminated")
 lines(X_3[,1],lm(X_3[,2]~X_3[,1])$fitted.value, pch=18, col="red", lty=1, lwd=1)
 abline(v=3, lty=2, col="blue")
 abline(2.0000000, 0.6666667, lty=3, col="black")
 legend(-5, 15, legend=c("LS",expression("T*"["RD"]),expression("T*"["PRD"])),col=c("red", "blue", "black"),
        lty=1:3, cex=0.7, title="Line types", text.font=4, bg='lightblue')
 text(11,1, "LS", col="red", cex=.8)
 text(4, 10, expression("T*"["RD"]), col="blue", cex=.8)
 text(13, 9, expression("T*"["PRD"]))