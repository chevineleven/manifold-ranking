function [D,sort,rank_score,query_output,positive_output,label_output]=manifold_ranking(distancematrix,query_input,positive_input,alpha,label_input,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: 
%distancematrix: distance matrix for manifold ranking.
%query_input: the index vector of query protein as positive samples.
%positive_input: the index vector of positive protein samples.
%label_input: the AC list of proteins.
%alpha: parameter in manifold ranking.
%flag: selction for manifold ranking algorithm:
%      0 - using converging form;
%      1 - using iteration form,iteration 100 times;


%Output:
%sort: a vector for the final order of proteins after manifold ranking.
%query_output: a index vector of query protein as positive samples after manifold ranking.
%positive_output: a index vector of positive protein samples after manifold ranking.
%label_output: the final AC list of proteins after manifold ranking.
%rank_score: the ranking score for each protein after manifold ranking.

%an example function call:[D,sort,rank_score,query_output,positive_output,label_output]=manifold_ranking(distancematrix,query_output,positive_output,0.5,label_output,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
query_length=length(query_input);
positive_length=length(positive_input);

[m,n]=size(distancematrix);

original_score=zeros(1,m);



for i=1:query_length

   original_score(query_input(i))=1; % set initial score
   
end;
original_score=original_score';

D=zeros(m,m);

for i=1:m
  for j=1:m

   D(i,i)=D(i,i)+distancematrix(i,j);
   
   end;
   D(i,i)=1/sqrt(D(i,i));
end;

L=D*distancematrix*D;

I=eye(m,m);

if (flag == 0)



rank_score=(1-alpha)*inv((I-alpha*L))*original_score;



elseif (flag == 1)



rank_score=original_score;
 for i=1:100

 
   rank_score=alpha*L*rank_score+(1-alpha)*original_score;

 end;



else return;


end;
f=rank_score;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

original_order=zeros(1,m);  %original_order is a vector to store the index of the original order of the protein in the input;
  for i=1:m
    original_order(i)=i;
    
  end;  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
original_query=zeros(1,m);
for i=1:query_length
   original_query(query_input(i))=1; %original_query is a 0/1 vector with 1 for query proteins and 0 for non-query proteins at each position
end;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

original_positive=zeros(1,m); %original_positive is a 0/1 vector with 1 for positive samples and 0 for negative samples at each position
for i=1:positive_length
   original_positive(positive_input(i))=1;
end;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:m 

 for j=(i+1):m
   
     if (f(j) > f(i))
     
         temp=f(i);
         f(i)=f(j);
         f(j)=temp;
         
         temp=original_order(i);
         original_order(i)=original_order(j);
         original_order(j)=temp;
         
         
         temp=original_query(i);
         original_query(i)=original_query(j);
         original_query(j)=temp;
         
         
         temp=original_positive(i);
         original_positive(i)=original_positive(j);
         original_positive(j)=temp;
         
     end;
end;

end;
rank_score=f;
sort=original_order;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate label_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label_output=cell(1,m);



for i=1:m
   label_output(i)=label_input(sort(i));
end;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Generate query_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
for i=1:m
   if (original_query(i)==1)
       count=count+1;
   end;
end;


query_output=zeros(1,count);


num=1;
for i=1:m
   if (original_query(i)==1)
       query_output(num)=i;
       num=num+1; 
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Generate positive_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
count=0;
for i=1:m
   if (original_positive(i)==1)
       count=count+1;
   end;
end;


positive_output=zeros(1,count);


num=1;
for i=1:m
   if (original_positive(i)==1)
       positive_output(num)=i;
       num=num+1; 
   end;
end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   



fid=fopen('MR_result.txt','w');
            
        
             for i=1:m
                 AC=char(label_output(i));
                 fprintf(fid,'%12s   ',AC);
                 
                 if (original_query(i)==1)
                   fprintf(fid,'-q   ');
                 end;  
                 
                 fprintf(fid,'%.10f',rank_score(i));
               
                 fprintf(fid,'\n');
             end;
             
fclose(fid);      
