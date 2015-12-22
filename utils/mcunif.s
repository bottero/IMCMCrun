N = 100000
a = 3000
b = 5000
mc = (a+b)/2
sig = 500
for (i in 1:N){
	u = rnorm(1,0,sig)
	if (mc[i]+u<a){
		mc = c(mc,mc[i]+u-a+b)
	}
	else if (mc[i] + u >b){
		mc = c(mc,mc[i]+u - b +a)
	}
	else {mc = c(mc,mc[i]+u)}
}
