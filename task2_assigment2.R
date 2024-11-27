b1=-1.2
b2=1.5
b3=1
b4=0.5
b5=1.3

a=exp(b1*2+b2*1+b3+4*b4)
b=exp(0)
c=exp(b1*2+b3*3+b5*1)
d=exp(1*b1+2*b4)
p1=a/(a+b+c+d)
p1

q=exp(2*b1+b2+b3+4*b4+b5)
w=exp(0)
e=exp(2*b1+b2+3*b4)
r=exp(b1+b4)
p2=q/(q+w+e+r)
p2

p=exp(b1+b2+b4)
o=exp(b1+2*b4+b5)
i=exp(3*b1+b2+b3+5*b4+b5)
u=exp(2*b1+b2+3*b4+b5)
p3=p/(p+o+i+u)
p3

lp=exp(2*b1+b2+4*b4)
lo=exp(3*b1+b2+b3+6*b4+b5)
li=exp(b1+2*b4)
lu=exp(b1+b2+2*b4)
p4=lp/(lp+lo+li+lu)
p4
