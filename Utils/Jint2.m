function Q = Jint2(x1,x2,y1,y2,xx,yy,f,fun,Qw)
x1= reshape(x1,1,1,[]);
x2= reshape(x2,1,1,[]);
y1= reshape(y1,1,1,[]);
y2= reshape(y2,1,1,[]);

xxm=bsxfun(@plus,x1,bsxfun(@times,xx,((x2-x1)./(2*f))));
yym=bsxfun(@plus,y1,bsxfun(@times,yy,((y2-y1)./(2*f))));
Q= mtimesx(Qw,mtimesx(fun(xxm,yym),Qw')).*(x2-x1).*(y2-y1)./4;
end

