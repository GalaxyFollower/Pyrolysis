function R=Gjinv(Q)
%Gjinv
%Matrix inverse using Gauss-Jordan elimination.
%Option of rational number format. 
%Option of all steps.
%Calling format: Gjinv(A)

%Can use given tolerance of 1e-14
%or change to own tolerance
%MATLAB computes to about 16 decimal digits

%Copyright Gareth Williams, Stetson University 
%gwilliam@stetson.edu, http://www.stetson.edu/~gwilliam
%Accompanies "Linear Algebra with Applications" by Gareth Williams

tol=1e-40;
[n,n]=size(Q);
P=[Q eye(n)];
[n,m]=size(P);
flag='T';

%find a pivot
if flag=='T'
 j=1;
 for i=1:n,
  if j <= m
   found=0;
   if abs(P(i, j)) <= tol  
     while (found == 0)
         
%search for a leading one and interchange rows if necessary
      for s=i:n,
       if (abs(P(s, j)) > tol) 
		  if  (found == 0)
        found=1;
         if s~=i
          for r=1:m,
           temp=P(i, r);
           P(i, r)=P(s, r);
           P(s, r) = temp;
          end

         end
        end
       end
      end

      if (found==0) 
       if (j <= n)
        flag='F'; %inverse does not exist
        found=1;
       end
      end

      if j>m  
       found=1;   % to exit while loop
      end 
     end %while

     if  j > m
      found = 0;
     end
    else
     found = 1;
   end  

%normalize leading element in row changing the rest of the row accordingly
   if flag=='T'
    if found == 1 
     k=i;
     if (P(k, j) ~= 1) 
      if (abs(P(k, j)) > tol)
       y = P(i, j);
       
        P(k, j:m) = P(k, j:m)/y ;
       

      end
     end
     for r=1:n,
      if (abs(P(r, j)) >tol) 
       if (r ~= i)
        z=P(r, j);
        P(r, j:m)=P(r, j:m) - z * P(i, j:m);
       end
      end
     end
    end
    j = j + 1;
   end
  end

 end %if flag 'T'

  R=P(:,n+1:end);

 if flag == 'F';
  disp('-Inverse Does Not Exist-')
 end 
end
