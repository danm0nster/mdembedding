% Run these tests with: runtests('mdDelayTest.m')

function tests = mdDelayTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add the directory above the test directory to search path
pathToAdd = strcat(pwd(), '/..');
addpath(pathToAdd, '-end');
% path
end

function teardownOnce(testCase)
disp("Teardown called. No action taken here.")
end

function testUniVariate(testCase)
x = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5]';
actTau = mdDelay(x, 'plottype', 'none');
expTau = 3;
verifyEqual(testCase, actTau, expTau);
end

function testBiVariate(testCase)
x = [1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5]';
y = [7, 6, 5, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 2, 3, 4, 5, 6]';
actTau = mdDelay([x,y], 'plottype', 'none');
expTau = 2.5;
verifyEqual(testCase, actTau, expTau);
end

function testMultivariateSine(testCase)
nDataPoints = 100;
intervalLength = 2*pi;
nVariables = 5;
x = zeros(nDataPoints, nVariables);
for var=1:nVariables
    x(:,var) = sin((1:nDataPoints)/nDataPoints*intervalLength + var);
end
actTau = mdDelay(x, 'plottype', 'none');
expTau = 2.8;
verifyEqual(testCase, actTau, expTau);
end