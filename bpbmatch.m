function [P, iter] = bpbmatch(W,b,draw);

%uses simplified belief propagation to compute optimal b-matching
%W is a weighted adjacency matrix. 
%
%[P, iter] = bpbmatch(W,b,draw);
%
%W = weight matrix
%b = degree of b-matching
%draw = 1 if graphical output is desired, 
%       0 if not (1 is slower but 
%         interesting to watch for small problems)
%P = output adjacency matrix
%iter = number of iterations to completion
%Copyright 2006 Bert Huang.
%Send comments and questions to bert at cs.columbia.edu

maxiter = 999;
oscmode=0;
count=0;

if (nargin==2)
    draw = 0;
end

P = zeros(size(W));

N = size(W,1);

W = exp(W); %exponentiate weight matrix

alpha = ones(N); %messages from alphas to betas
beta = ones(N); %message from betas to alphas

foundanswer = logical(0);

workvec = zeros(N,1);
inds = zeros(N,1);

iter=1;

prevstates = [];

code = rand(N);
basestate=0;

while (foundanswer<=0)
    %iterate message passing
    
    %update alpha
    for i=1:N
        %compute work vector
        workvec = W(:,i).*alpha(:,i);

        inds = 1:b+1;
        [junk, minind] = min(workvec(inds));
        for j=b+2:N
            if (workvec(j)>workvec(inds(minind)))
                inds(minind) = j;
                [junk, minind] = min(workvec(inds));
            end
        end
        auxvec = workvec(inds);
        auxvec(minind) = inf;
        [junk, truebth] = min(auxvec);
       
        for j=1:N
            if (workvec(j)>=workvec(inds(truebth))) %check if j is in the top b
                beta(i,j) = W(j,i)/workvec(inds(minind));
            else
                beta(i,j) = W(j,i)/workvec(inds(truebth));
            end
        end
    end
       
    for i=1:N
        %compute work vector
        workvec = W(i,:)'.*beta(:,i);

        inds = 1:b+1;
        [junk, minind] = min(workvec(inds));
        for j=b+2:N
            if (workvec(j)>workvec(inds(minind)))
                inds(minind) = j;
                [junk, minind] = min(workvec(inds));
            end
        end
        auxvec = workvec(inds);
        auxvec(minind) = inf;
        [junk, truebth] = min(auxvec);
       
        for j=1:N
            if (workvec(j)>=workvec(inds(truebth))) %check if j is in the top b
                alpha(i,j) = W(i,j)/workvec(inds(minind));
            else
                alpha(i,j) = W(i,j)/workvec(inds(truebth));
            end
        end
    end
    
    %compute "beliefs"

    P = zeros(N);
    for i=1:N
        workvec = W(i,:)'.*beta(:,i);

        inds = 1:b;
        [junk, minind] = min(workvec(inds));
        for j=b+1:N
            if (workvec(j)>workvec(inds(minind)))
                inds(minind) = j;
                [junk, minind] = min(workvec(inds));
            end
        end
        
        P(i*ones(1,b), inds)=1;
    end
    if (draw==1)
        pvec = sum(P);
        pgreen = repmat(pvec==b, N,1);
        pvec = sum(P');
        pred = repmat(pvec'==b,1,N);
        Pcolor=zeros(N,N,3);
        Pcolor(:,:,1) = pred.*P;
        Pcolor(:,:,2) = pgreen.*P;
        Pcolor(:,:,3) = P;

        subplot(311);
        imagesc(W); title('W');
        subplot(312);
        imagesc(log(1+alpha)); title('mu'); xlabel(num2str(iter));
        subplot(313);
        imagesc(Pcolor);
        colormap gray;
        drawnow;
    end

    %check for answer
    if (iter>N && sum(abs(sum(P)-b*ones(1,N))+abs(sum(P')-b*ones(1,N)))==0)
        foundanswer=foundanswer+1;
    else
        foundanswer=foundanswer-.1;
    end

    if (iter>=maxiter)
        disp('*****************************');
        disp('WARNING: Reached maximum iterations');
        disp('*****************************');
        foundanswer=1;
    end
        
    iter = iter+1;
end
