function [ output_args ] = createKnowledgeLib( typeClass, num_attr )
%% CREATEKNOWLEDGELIB: the function to create the knowledge library
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  typeClass  --> the type of species
         num_attr   --> the number of attributes
 Output: output_args --> the universal set
 [Tips]: "typeClass" the type of species({'setosa'}; {'versicolor'}; {'virginica'})
         "num_attr"  the number of attributes(for IRIS dataset is 4)
         (SL, SW, PL, PW)
%}

output_args = []; % allocate the variant

num_class = length(typeClass);
tmp_fe = cell(1, num_attr);

for n1 = 1:num_class % for SL
    tmp_fe(1) = typeClass(n1);
    for n2 = 1:num_class % for SW
        tmp_fe(2) = typeClass(n2);
        for n3 = 1:num_class % for PL
            tmp_fe(3) = typeClass(n3);
            for n4 = 1:num_class % for PW
                tmp_fe(4) = typeClass(n4);
                output_args = [output_args; {tmp_fe}];
            end
        end
    end
end
% the number of planes    
% for n = 1:length(lib.num_planes_select)
%     % the all possible type selection
%     for m = 1:length(lib.type_select)
%         % the all possible curvature selection
%         for l = 1:length(lib.curvature_select)
%             % setting 
%             num_planes = lib.num_planes_select(n);
%             type = lib.type_select(m);
%             curvature = lib.curvature_select(l);
%             
%             tmp_element = createBasicElement(num_planes, type, curvature);
%             output_args = [output_args, tmp_element];
%         end
%     end
% end

end

