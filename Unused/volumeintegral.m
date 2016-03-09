function [volint]=volumeintegral(X,Y,upwind)
[a]=polynomiame(X,Y,2);
a=[a zeros([size(a,1) 2])];
volint=integratethem(a,X,upwind);




