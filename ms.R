library(ggplot2)

# dinvgamma: density function of the inverse-gamma distribution
dinvgamma = function(x, shape = 1, rate = 1, scale = 1/rate, log = FALSE) {
    # return( rate^shape / gamma(shape) * x^(-shape-1) * exp(-rate/x) )
    logval = shape*log(rate) - lgamma(shape) - (shape+1)*log(x) - rate/x
    if (log)
	return(logval)
    else
	return(exp(logval))
}

# pinvgamma: cumulative distribution function of the inverse-gamma distribution
pinvgamma = function(q, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( pgamma(1/q,shape,rate,scale,!lower.tail,log.p) )
}

# qinvgamma: quantile function of the inverse-gamma distribution
qinvgamma = function(p, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( 1 / qgamma(p,shape,rate,scale,!lower.tail,log.p) )
}

# rinvgamma: sample a vector of inverse-gamma random variables
rinvgamma = function(n, shape = 1, rate = 1, scale = 1/rate) {
    return( 1 / rgamma(n,shape,rate,scale) )
}

# base graphics
baseplot <- function(...){
  plot ( x=c(0,1.3), y=c(0,1.3), type="n", bty="l",
         xlab = expression(sigma[e]^2), ylab = expression(sigma[s]^2),
         line=2,
         xaxt="n", yaxt="n",
         ...
       )
}

line1 <- function ( label=TRUE, ... ){
  lines ( x=c(0,.4), y=c(.7,0), ... )
  if(label) text ( x=.23, y=.3, expression(j==1) )
}
line2 <- function ( label=TRUE, ... ){
  lines ( x=c(0,.1), y=c(1,0), ... )
  if(label) text ( x=.1, y=.1, expression(j==2) )
}
line3 <- function ( label=TRUE, ... ){
  lines ( x=c(0,1), y=c(.6,0), ... )
  if(label) text ( x=.5, y=.3, expression(j==3) )
}
line0 <- function ( label=TRUE, ... ){
  lines ( x=c(.6,.6), y=c(1,0), ... )
  if(label) text ( x=.6, y=.7, expression(j==0) )
}
linev <- function ( label=TRUE, ... ){
  lines ( x=c(.6,.6), y=c(1,0), ... )
  if(label) text ( x=.61, y=.4, expression(j==v) )
}
lineh <- function ( label=TRUE, ... ){
  lines ( x=c(0,1), y=c(.6,.6), ... )
  if(label) text ( x=.3, y=.62, expression(j==h) )
}


boxA <- function(){
  lines ( x=c(.1, .15, .15, .1, .1), y=c(.2, .2, .25, .25, .2) )
  text ( x=.125, y=.225, "A", cex=1.5)
}
boxB <- function(){
  lines ( x=c(.2, .25, .25, .2, .2), y=c(.5, .5, .55, .55, .5) )
  text ( x=.225, y=.525, "B", cex=1.5)
}
boxC <- function(){
  lines ( x=c(.3,.35,.35,.3,.3), y=c(.1,.1,.15,.15,.1) )
  text ( x=.325, y=.125, "C", cex=1.5)
}
boxb <- function(){
  lines ( x=c(.3,.35,.35,.3,.3), y=c(.1,.1,.15,.15,.1) )
  text ( x=.325, y=.125, expression(italic(b)), cex=1.5)
}

# ggplot2 version
df <- data.frame ( x=c(0,1.3), y=c(0,1.3) )
p <- ggplot ( df, aes(x=x,y=y) ) +
  scale_x_continuous ( limits = c(0,1.3),
                       name = expression(sigma[e]^2),
                       breaks = NULL
                     ) +
  scale_y_continuous ( limits = c(0,1.3),
                       name = expression(sigma[s]^2),
                       breaks = NULL
                     ) +
  theme_bw()
x1 <- c(0,.4); x1ann <- 0.23; y1 <- c(.7,0); y1ann <- 0.3
x2 <- c(0,.1); x2ann <- 0.1; y2 <- c(1,0); y2ann <- 0.1
x3 <- c(0,1); x3ann <- 0.5; y3 <- c(.6,0); y3ann <- 0.3
xv <- c(.6,.6); xvann <- 0.61; yv <- c(0,1); yvann <- 0.4
xh <- c(0,1); xhann <- 0.3; yh <- c(.6,.6); yhann <- 0.62
B1 <- data.frame ( x=c(0,1,1,0,0), y=c(0,0,1,1,0) )
B1annx <- c (.2,.97); B1anny <- c(.97,.2)
# B1small <- data.frame ( x=c(1,.6,.6,0,0,2/47,.1,1),
#                         y=c(0,.24,1,1,.6,27/47,0,0)
#                       )
B1small <- data.frame ( x=c(1,.6,.6,.4/7,4/110,0,0,.4/9.4,.1,1),
                        y=c(0,.24,.6,.6,7/11,1,.6,5.4/9.4,0,0)
                      )
B <- data.frame ( x=c(0,1.2,1.2,0,0), y=c(1.2,1.2,0,0,1.2) )
Bannx <- c (.2,1.17); Banny <- c(1.17,.2)

# Figure "oneline"
pdf("figs/oneline.pdf",width=5,height=5)
# baseplot()
# line1()
# text ( x=0.1, y=0.2, "+", cex=2 )
# text ( x=0.3, y=0.4, "-", cex=2 )
p + annotate ( "line", x=x1, y=y1 ) +
  annotate ( "text", x=0.1, y=0.2, size=12, label="+") +
  annotate ( "text", x=0.3, y=0.4, size=12, label="-")
dev.off()

# Figure "multilines"
pdf("figs/multilines.pdf")
baseplot()
line1()
line2()
line3()
line0()
boxA()
boxB()
boxC()
p + l1 + l2 + l3 + lv + lh
dev.off()

#Figure "boundingbox"
pdf("figs/boundingbox.pdf")
# baseplot()
# line1(lty=3)
# line2(lty=3)
# line3(lty=3)
# #line0(lty=3)
# linev(lty=3)
# lineh(lty=3)
# lines ( x=c(1,.6,.6,0,0,2/47,.1,1), y=c(0,.24,1,1,.6,27/47,0,0), lty=2 )
# lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=2 )
# text(x=.2, y=.98, expression(B[1]))
# text(x=.98, y=.2, expression(B[1]))
# points ( x=0.5, y=1.2, pch=19 ); text ( x=.55, y=1.2, expression(bolditalic(p)[1]))
# points ( x=0.5, y=1.0, pch=19 ); text ( x=.55, y=1.05, expression(bolditalic(p)[1]^"*"))
# points ( x=1.05, y=0.5, pch=19 ); text ( x=1.1, y=0.5, expression(bolditalic(p)[2]))
# points ( x=1.0, y=0.5, pch=19 ); text ( x=0.95, y=0.5, expression(bolditalic(p)[2]^"*"))
# points ( x=1.05, y=1.1, pch=19 ); text ( x=1.05, y=1.15, expression(bolditalic(p)[3]))
# points ( x=1, y=1, pch=19 ); text ( x=1.05, y=1, expression(bolditalic(p)[3]^"*"))
p + annotate ( "line", x=x1, y=y1, lty=3 ) +
  annotate ( "line", x=x2, y=y2, lty=3 ) +
  annotate ( "line", x=x3, y=y3, lty=3 ) +
  annotate ( "line", x=xv, y=yv, lty=3 ) +
  annotate ( "line", x=xh, y=yh, lty=3 ) +
  annotate ( "text", x=x1ann, y=y1ann, size=4, parse=TRUE, label="j==1" ) +
  annotate ( "text", x=x2ann, y=y2ann, size=4, parse=TRUE, label="j==2" ) +
  annotate ( "text", x=x3ann, y=y3ann, size=4, parse=TRUE, label="j==3" ) +
  annotate ( "text", x=xvann, y=yvann, size=4, parse=TRUE, label="j==s[Z]+1" ) +
  annotate ( "text", x=xhann, y=yhann, size=4, parse=TRUE, label="j==s[Z]+2" ) +
  geom_path ( data=B1, size=.7 ) +
  annotate ( "text", x=B1annx, y=B1anny, size=4, parse=TRUE, label="B[1]") +
  annotate ( "point", x=c(.5,.5, 1.05, 1, 1.05, 1), y=c(1.2,1,.5,.5,1.1,1) ) +
  annotate ( "text", x=c(.55,.55,1.1,.95,1.05,1.05),
                     y=c(1.2,1.05,.5,.5,1.15,1), parse=TRUE,
             label = c ( "p[1]", "widetilde(p[1])",
                         "p[2]", "widetilde(p[2])",
                         "p[3]", "widetilde(p[3])" )
           )
dev.off()

#Figure "smallboundingregion"
pdf("figs/smallboundingregion.pdf")
# baseplot()
# line1(lty=3)
# line2(lty=3)
# line3(lty=3)
# linev(lty=3)
# lineh(lty=3)
#lines ( x=c(1,.6,.6,0,0,2/47,.1,1), y=c(0,.24,1,1,.6,27/47,0,0), lwd=2 )
# lines ( x=c(1,1,.6,.6,0,0,2/47,.1,1), y=c(0,.6,.6,1,1,.6,27/47,0,0), lwd=2 )
# lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lty=2 )
p + annotate ( "line", x=x1, y=y1, lty=3 ) +
  annotate ( "line", x=x2, y=y2, lty=3 ) +
  annotate ( "line", x=x3, y=y3, lty=3 ) +
  annotate ( "line", x=xv, y=yv, lty=3 ) +
  annotate ( "line", x=xh, y=yh, lty=3 ) +
  annotate ( "text", x=x1ann, y=y1ann, size=4, parse=TRUE, label="j==1" ) +
  annotate ( "text", x=x2ann, y=y2ann, size=4, parse=TRUE, label="j==2" ) +
  annotate ( "text", x=x3ann, y=y3ann, size=4, parse=TRUE, label="j==3" ) +
  annotate ( "text", x=xvann, y=yvann, size=4, parse=TRUE, label="j==s[Z]+1" ) +
  annotate ( "text", x=xhann, y=yhann, size=4, parse=TRUE, label="j==s[Z]+2" ) +
  geom_path ( data=B1, lty=3 ) +
  geom_path ( data=B1small, size=.7 )
dev.off()


#Figure "boundingbox2"
pdf("figs/boundingbox2.pdf")
# baseplot()
# line1(lty=3,label=FALSE)
# line2(lty=3,label=FALSE)
# line3(lty=3,label=FALSE)
# line0(lty=3,label=FALSE)
# lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=1 )
# text (x=.97, y=.5, expression(B[1]))
# text (x=.3, y=.97, expression(B[1]))
# points ( x=0.25, y=.7, pch=19 ); text ( x=.25, y=.75, expression(bolditalic(p)[4]))
# lines (x=c(0,1.2,1.2,0,0), y=c(1.2,1.2,0,0,1.2), lwd=2)
# text (x=1.17, y=.5, expression(B))
# text (x=.3, y=1.17, expression(B))
# axis(side=1,at=1.2,labels=expression(sigma[e]^"2*"))
# axis(side=2,at=1.2,labels=expression(sigma[s]^"2*"))
# lines(x=c(1.2,1.2),y=c(1.2,1.3),lty=2)
# lines ( x=c(.6,.6), y=c(1,1.3), lty=2)
# lines(x=c(0,0),y=c(1.2,1.3),lty=2)
# points(x=1.25,y=1.25,pch=19); text(x=1.25,y=1.25,expression(bolditalic(q)[1]),pos=4)
# points(x=1.2,y=0,pch=19); text(x=1.2,y=0,expression(bolditalic(q)[1]^"*"),pos=4)
#points(x=1,y=0,pch=19); text(x=1,y=0.05,expression(bolditalic(q)[1]^"**"),pos=4)
# points (x=.8, y=1.25, pch=19); text (x=.8, y=1.25, expression(bolditalic(q)[2]),pos=4)
# points (x=.6, y=1.2, pch=19); text(x=.6,y=1.225, expression(bolditalic(q)[2]^"*"),pos=4)
#points (x=.6, y=1, pch=19); text(x=.6,y=1.025, expression(bolditalic(q)[2]^"**"),pos=4)
# points(x=.2,y=1.25,pch=19); text(x=.2, y=1.25, expression(bolditalic(q)[3]),pos=4)
# points(x=.0,y=1.2,pch=19); text(x=0, y=1.225, expression(bolditalic(q)[3]^"*"),pos=4)
#points(x=.0,y=1,pch=19); text(x=0, y=1.025, expression(bolditalic(q)[3]^"**"),pos=4)
p + annotate ( "line", x=x1, y=y1, lty=3 ) +
  annotate ( "line", x=x2, y=y2, lty=3 ) +
  annotate ( "line", x=x3, y=y3, lty=3 ) +
  annotate ( "line", x=xv, y=c(0,1.2), lty=3 ) +
  annotate ( "line", x=c(0,1.2), y=yh, lty=3 ) +
  annotate ( "text", x=x1ann, y=y1ann, size=4, parse=TRUE, label="j==1" ) +
  annotate ( "text", x=x2ann, y=y2ann, size=4, parse=TRUE, label="j==2" ) +
  annotate ( "text", x=x3ann, y=y3ann, size=4, parse=TRUE, label="j==3" ) +
  annotate ( "text", x=xvann, y=yvann, size=4, parse=TRUE, label="j==s[Z]+1" ) +
  annotate ( "text", x=xhann, y=yhann, size=4, parse=TRUE, label="j==s[Z]+2" ) +
  geom_path ( data=B1, lty=2 ) +
  annotate ( "text", x=B1annx, y=B1anny, size=4, parse=TRUE, label="B[1]") +
  annotate ( "text", x=Bannx, y=Banny, size=4, parse=TRUE, label="B") +
  geom_path ( data=B, size=.7 ) +
  annotate ( "point",
             x=c(.2, .8, 1.25, 1.25, 0, .6, 1.2, 1.2, .3),
             y=c(1.25, 1.25, 1.25, .2, 1.2, 1.2, 0, .6, .8)
           ) +
  annotate ( "text",
             x=c(.2, .8, 1.28, 1.28, .02, .6, 1.24, 1.24, .33),
             y=c(1.28, 1.28, 1.28, .2, 1.24, 1.24, .6, 0, .8),
             parse=TRUE,
             label = c ( "q[1]", "q[2]", "q[3]", "q[4]",
                         "widetilde(q[1])", "widetilde(q[2])",
                         "widetilde(q[3])", "widetilde(q[4])", "p"
                       )
           ) +
  scale_x_continuous ( name = expression(sigma[e]^2),
                       breaks=1.2,
                       labels=expression(widetilde(sigma[e]^2))
                     ) +
  scale_y_continuous ( name = expression(sigma[s]^2),
                       breaks=1.2,
                       labels=expression(widetilde(sigma[s]^2))
                     )
dev.off()

# Figure "abcboxes"
pdf("figs/abcboxes.pdf")
baseplot()
line1()
boxA()
boxB()
boxC()
dev.off()

# Figure "boxb"
pdf("figs/boxb.pdf")
baseplot()
lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=2 )
text(x=.2, y=.98, expression(B[1]))
text(x=.98, y=.2, expression(B[1]))
line1(lty=3)
line2(lty=3)
line3(lty=3)
boxb()
dev.off()


#Figure "contour box"
pdf("figs/contourbox.pdf")
baseplot(usr=c(0,1.2,0,1.2))
line1(label=F,lty=2)
line2(label=F,lty=2)
line3(label=F,lty=2)
line0(label=F,lty=2)
arrows ( 1, 0, 1, 1.04 ); text ( x=.98, y=.5, expression(bolditalic(l)[1]))
points ( x=1, y=0, pch=19 ); text ( x=.97, y=0, expression(bolditalic(p)[1]))
arrows ( 0, 1, 1, 1 ); text ( x=.5, y=.98, expression(bolditalic(l)[2]))
points ( x=0, y=1, pch=19 ); text ( x=.022, y=.98, expression(bolditalic(p)[2]))
points ( .6, 1, pch=19 ); text ( x=.62, y=.98, expression(bolditalic(p)[3]))
dev.off()
