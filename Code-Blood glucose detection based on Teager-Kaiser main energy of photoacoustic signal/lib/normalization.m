function [ output_args, max_v, min_v ] = normalization( input_args )
%NORMALIZATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%{
   input_args: 1D vector
%}
max_v = max(input_args);
min_v = min(input_args);

output_args = (input_args - min_v)/(max_v - min_v);

end

