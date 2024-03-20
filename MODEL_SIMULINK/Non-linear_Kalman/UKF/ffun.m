function Xnew=ffun(X) % 状态方程函数 
Xnew(1,1)=X(1)+1; 
Xnew(2,1)=X(2)+sin(0.1*X(1));