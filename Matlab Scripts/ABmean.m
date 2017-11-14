%Calculates Alpha, Beta from mean of raw distance data
%Called as ABmean(min,max,step)

% Copyright 2017 Advanced Detector Centre, Department of Electronic and
% Electrical Engineering, University of Sheffield, UK. 
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function f=ABmean(min,max,step)
n=1+(max-min)/step;

fileID=fopen('ab.txt','w');
for k=1:n
    bias=min+step*(k-1);
    efile=sprintf('%gepdf.txt',bias);
    hfile=sprintf('%ghpdf.txt',bias);
    edata=dlmread(efile);
    hdata=dlmread(hfile);
    alpha=1/mean(edata(:,2));
    beta=1/mean(hdata(:,2));
    fprintf(fileID,'%g %g %g\r\n', bias,alpha,beta);
end
fclose(fileID);
end
    