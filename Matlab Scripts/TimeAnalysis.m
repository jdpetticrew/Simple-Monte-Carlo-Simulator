% Returns Mean time to breakdown and Jitter for a SMC dataset
% [Meantime,Jitter]

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
function t=TimeAnalysis(time)
time=time/1e-12;
meantime=mean(time);
jitter=FWHM2(time);
t=[meantime,jitter];
end

function f= FWHM2(time)
h=histogram(time);
h.BinWidth=0.1;
[M,I]=max(h.Values(:));
TimeShift=abs(h.Values-M/2);
[~,K]=min(TimeShift(1:I));
[~,Y]=min(TimeShift(I+1:numel(TimeShift)));
f=((I+Y)-K)*h.BinWidth;
end