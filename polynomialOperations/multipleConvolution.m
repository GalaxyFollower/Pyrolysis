function C=multipleConvolution(A, B)
    if size(A,1)~=size(B,1)
        error(' A and B must have the same length');
    end
    C=zeros([size(A,1) size(A,2)+size(B,2)-1]);
    for i=1:size(A,1)
        C(i,:)=conv(A(i,:),B(i,:));
    end
end