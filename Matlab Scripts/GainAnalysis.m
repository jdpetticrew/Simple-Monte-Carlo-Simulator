% Calculates Noise and Gain from SMC Gain Out Data [Gain, Noise]

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

function out = GainAnalysis(Gain)
Gain=Gain/1.6e-19;
MeanGain=mean(Gain);
Gain2=Gain.*Gain;
MeanGain2=mean(Gain2);
Mean2=MeanGain*MeanGain;
f=MeanGain2/Mean2;
m=MeanGain;
out=[m,f];
end
