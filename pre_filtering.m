function [sort,vector_matrix,distancematrix,query_output,positive_output,label_output]=pre_filtering(filename,sigma,query_input,throwaway,positive_input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: 
%filename: the name of input protein dataset file.
%sigma: parameter in manifold ranking to calculate the weight between point i and j.
%query_input: the index vector of query protein as positive samples.
%throwaway: the number of proteins to be fiterted out.
%positive_input: the index vector of positive protein samples.

%Output:
%sort: a vector for the final order of proteins according to their weight to nearest query.
%query_output: a index vector of query protein as positive samples after filtering.
%positive_output: a index vector of positive protein samples after filtering.
%label_output: the final AC list of proteins according to their weight to nearest query.
%vector_matrix: vector matrix to store the protein samples in weight order after filtering.
%distancematrix: distance matrix for manifold ranking.

%an example function call:[sort,vector_matrix,distancematrix,query_output,positive_output,label_output]=pre_filtering('ex.sub.dat',1.25,[1:52],0,[1:52]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load protein dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_matrix=importdata(filename);
sub_matrix=total_matrix.data;
[m,feature_length]=size(sub_matrix);
sub_matrix(:,feature_length)=[];
AC_label=total_matrix.textdata;  % AC_label is a cell vector to store the accession number of proteins
AC_label=AC_label';
original_vector_matrix=sub_matrix;

label_length=length(AC_label);
query_length=length(query_input);
positive_length=length(positive_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

original_order=zeros(1,m);  %original_order is a vector to store the index of the original order of the protein in the input;
  for i=1:m
    original_order(i)=i;
  end;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

original_query=zeros(1,m);  %original_query is a 0/1 vector with 1 for query proteins and 0 for non-query proteins at each position
for i=1:query_length
   original_query(query_input(i))=1;
end;   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  

original_positive=zeros(1,m); %original_positive is a 0/1 vector with 1 for positive samples and 0 for negative samples at each position
for i=1:positive_length
   original_positive(positive_input(i))=1;
end;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% nearest neighbor rule for weight assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


weight_vector=zeros(1,m); %weight_vector contains the weight for each sample to the nearest query. It will be used for further sorting
for i=1:m
    minweight=0;
    weight=100;
    vector1=sub_matrix(i,:);
    
    for j=1:query_length
       
       vector2=sub_matrix(query_input(j),:);
       
          
          tempweight=exp(-((norm(vector1-vector2))^2)/(2*sigma^2));
          
          
          
          if (tempweight>minweight)
              weight=tempweight;
              minweight=tempweight;
          end;
          
          if (i==query_input(j))
            weight=100;  % if the point is a query sample, its weight is always the max
            break;
          end;
          
              
    end;
     
    weight_vector(i)=weight;
    
end;            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:m 

 for j=(i+1):m
   
     if (weight_vector(j) > weight_vector(i))
     
         temp=weight_vector(i);
         weight_vector(i)=weight_vector(j);
         weight_vector(j)=temp;
         
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
sort=original_order; % sort is a vector for the final order of proteins according to their weight to nearest query
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate label_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

label_output=cell(1,m);

for i=1:m
   label_output(i)=AC_label(sort(i)); % label_output is the final AC list of proteins according to their weight to nearest query
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate vector_matrix to store the protein dataset after filtering in order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vector_matrix=zeros(m-throwaway,feature_length-1);

for i=1:(m-throwaway)

    vector_matrix(i,:)=original_vector_matrix(sort(i),:);
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate distancematrix for manifold ranking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


distancematrix=zeros((m-throwaway),(m-throwaway));

      
for i=1:(m-throwaway)
   
        vector1=vector_matrix(i,:);
        
       for j=(i+1):(m-throwaway) 
      
           vector2=vector_matrix(j,:);
           
           distancematrix(i,j)=exp(-((norm(vector1-vector2))^2)/(2*sigma^2));
           distancematrix(j,i)=distancematrix(i,j);
       end;
end;           





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('filtering_result.txt','w');
            
        
             for i=1:m
                 AC=char(label_output(i));
                 fprintf(fid,'%12s ',AC);
                 
                 if (original_query(i)==1)
                   fprintf(fid,'-q');
                 end;  
                 
               
                 fprintf(fid,'\n');
             end;
             
fclose(fid);                