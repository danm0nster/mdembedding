% Run these tests with: runtests('autoMITest.m')

function tests = autoMITest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add the directory above the test directory to search path
pathToAdd = strcat(pwd(), '/..');
addpath(pathToAdd, '-end');
% path
end

function testReference(testCase)
x = sin(1:100)/4/pi';
res = autoMI(x, 10, 10);
reference = mi(x, 10, 10, 'silent');
ref = squeeze(reference(1,1,:));
verifyEqual(testCase, res, ref, 'AbsTol', 0.1);
end