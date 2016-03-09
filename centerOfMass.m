function cm=centerOfMass(a,I,r)
[np m]=size(a);
q=m-1;

a_den=[a./repmat(q+4:q+1,[np 1]) zeros([np q+1 ])] ;
a_num=[a./repmat(q+3:q,[np 1]) zeros([np q])];



for i=1:length(I)-1
	j=find(r>=I(i) & r<=I(i+1));
    n=length(j);
    cm(j(1:n-1))=(polyval(a_den(i,:),r(j(2:n))) ...
        -polyval(a_den(i,:),r(j(1:n-1))))./(...
        polyval(a_num(i,:),r(j(2:n)))-polyval(a_num(i,:),r(j(1:n-1))));    
end

end




