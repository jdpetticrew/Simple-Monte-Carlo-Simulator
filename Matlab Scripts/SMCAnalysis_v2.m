% Automation of TimeAnalysis and GainAnalysis
%Called as SMCAnalysis_v2() reads Simple Monte Carlo Simulator files for 
%it's variables.

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


function f=SMCAnalysis_v2()
bias=dlmread('bias_input.txt');
[numvoltages,~]=size(bias);
num=cast(numvoltages,'int16');
results= zeros(num,5);

for k=1:numvoltages
    voltage=bias(k);
	results(k,1)=voltage;
	try
		% Trys to analyse gain data from SMC
		GainFile = sprintf('%ggain_out.txt',voltage);
		GainData= dlmread(GainFile);    
		pass1=GainData(:,2);
		dump=GainAnalysis(pass1);
		results(k,2)=dump(1);
		results(k,3)=dump(2);
	catch err
		% Displays missing file
		GainFile = sprintf('%ggain_out.txt does not exist',voltage);
		disp(GainFile)
	end
	
	try
		%Trys to analyse noise data from SMC.
		TimeFile = sprintf('%gtime_to_breakdown.txt',voltage);
		TimeData = dlmread(TimeFile);
		dump=TimeAnalysis(TimeData(:,2));
		results(k,4)=dump(1);
		results(k,5)=dump(2);
	catch err
		% Displays missing file
		TimeFile = sprintf('%gtime_to_breakdown.txt does not exist',voltage);
		disp(TimeFile)
	end
end

%Output results to text file.
fileID=fopen('results_matlab.txt','w');
fprintf(fileID,'%3s %5s %6s %9s %7s\r\n','V','Gain','Noise','MeanTime','Jitter');

for j=1:numvoltages
    fprintf(fileID,'%f %f %f %f %f\r\n',results(j,1),results(j,2),results(j,3),results(j,4),results(j,5));
end
fclose(fileID);

%Returns Results
f=results;
end