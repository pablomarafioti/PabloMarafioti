pan=function(x1,x2,alt,B=10000)
{pan.stat=function(x1,x2)
 {scarto.mediana=function(data)
  {a=sort(abs(data-median(data)))
  if (length(a)%%2)
  a[a!=0 | duplicated(a)]
  else a}
  panstat=(log(mean(scarto.mediana(x1)))-log(mean(scarto.mediana(x2))))/(var(scarto.mediana(x1))/length(scarto.mediana(x2))/(mean(scarto.mediana(x1)))^2+var(scarto.mediana(x2))/length(scarto.mediana(x1))/(mean(scarto.mediana(x2)))^2)^0.5
 }

n1=length(x1)
n2=length(x2)
x1=x1-mean(x1)
x2=x2-mean(x2)
pan.obs=pan.stat(x1,x2)
pan.perm=vector(,B)
x=c(x1,x2)

for (b in 1:B)
{x.perm=sample(x)
 x1.perm=x.perm[1:n1]
 x2.perm=x.perm[(n1+1):(n1+n2)]
 pan.perm[b]=pan.stat(x1.perm,x2.perm)}

if (alt=="greater") pvalue=length(pan.perm[pan.perm>=pan.obs])/B
if (alt=="less") pvalue=length(pan.perm[pan.perm<=pan.obs])/B
if (alt=="two.sided") pvalue=length(abs(pan.perm)[abs(pan.perm)>=abs(pan.obs)])/B

output=list(obs.value=pan.obs,p.value=pvalue)
}
