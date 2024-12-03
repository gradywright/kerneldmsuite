function varargout = kdmtest(varargin)
%KDMTEST   Run Kernel DM test suite.
%   KDMTEST executes all of the m-files found in the top level folders of the
%   directory $kerenldmsuiteroot/tests/. These m-files should return a scalar,
%   vector, or matrix of logical values. A test is deemed to pass if all the
%   returned values are logical true.  There is no functional output from
%   KDMTEST, but the data is piped to the command window with fprintf. Note
%   that KDMTEST will automatically close any open figure windows.
%
%   KDMTEST will time each of the tests, and print the time for each test
%   that passes. If any of the values returned by a test m-file are logical
%   false (or zero) then the test is deemed to have 'failed'. If a test
%   results in an error being thrown it is reported as 'crashed'. 
%
%   If all test from a particular directory pass, the total execution time is
%   returned.  If not, the number of failed tests is reported.
%
%   Similarly, if all tests pass in all directories the total execution time is
%   returned at the end of the test. If not, then a list of the failed tests 
%   is reported:
%
%   RESULTS = KDMTEST returns an N-by-3 cell-array with one row per test and 
%   columns:
%     (1) directory name, a string.     Example: 'sphere'
%     (2) test name, a string.          Example: 'test_rbffdlap'
%     (3) duration OR error flag. A value >0 indicates success and is equal to
%         the execution time of the test in seconds, while -1 == failure and
%         -2 == crash.

% This is a modified version of CHEBTEST from the chebfun package: 
% http://www.chebfun.org/

% Find directory in which kerneldmsuite was installed:
installDir = kdmroot();
quietMode = false;

% Set path to the tests/ subdirectory:
testsDir = fullfile(installDir, 'tests');

% Run all the test we can find in the tests/ folder.
testDirNames = dir(testsDir);       % List test directories.
testDirNames = {testDirNames.name}; % Keep only the directory names.
testDirNames = testDirNames(3:end); % Strip '.' and '..' directories.

% If testDirNames is empty here, something is wrong.
if ( isempty(testDirNames) )
    error('KDMSUITE:kdmtest:testsNotFound', ...
        ['Could not locate test directories. ' ...
         'Please check that Kernel DM Suite has been installed correctly.']);
end

% Store the current directory: (We will return here when we're done.)
currDir = pwd();

% Switch to the tests/ subdirectory:
cd(testsDir);

% Initialise storage for which directories pass and the time they take.
numDirs = length(testDirNames);
allResults = cell(0, 3);

% Loop over the test directories and run the tests in each one.
for k = 1:numDirs
    testDir = testDirNames{k};
    if ( exist(testDir, 'dir') )
        fprintf(['Running tests in ' testDir ':\n']);
        nextResults = runTestsInDirectory(testDir, quietMode);
        allResults = [allResults ; nextResults]; %#ok<AGROW>
    else
        warning('KDMSUITE:kdmtest:dirNotFound', ...
            'Test directory ''%s'' not found. Skipping.', testDir);
    end
    fprintf('\n');
end

durations = [allResults{:,3}];
if ( all(durations > 0) )
    % All tests have passed. Yay!
    
    % Note that the timing information is superfluous if only testing one dir.
    if ( numDirs > 1 )
        fprintf('All tests passed in %4.4fs\n', sum(durations));
    end

else
    % There's been a failure. :(
    
    % Note that we don't display this in quiet mode, as the failed tests are
    % already easily seen.
    if ( ~quietMode )
        indx = find(durations < 0);
        fprintf('The following tests failed or crashed:\n');
        for k = indx(:)'
            loc = fullfile(testsDir, allResults{k,1});
            fprintf('   %s\n', printTestFilename(allResults{k,1:2}, loc))
        end
    end

end

% Return to current directory.
cd(currDir);

% Return.
if ( nargout > 0 )
    varargout{1} = allResults;
end

end


function testResults = runTestsInDirectory(testDir, quietMode)
%RUNTESTSINDIRECTORY   Run all the tests in the given directory.
%   RUNTESTSINDIRECTORY(TESTDIR) will change the current working directory to
%   TESTDIR, locate all the *.m files within it using DIR, and execute each of
%   these in turn (in alphabetical order). If any of the tests crash (i.e.,
%   throw an error), it will be caught in a try-catch statement, and the error
%   will not be rethrown.
%
%   If all the tests in TESTDIR pass, then the total execution time for this
%   directory is also printed to screen.
%
%   RUNTESTSINDIRECTORY returns a cell-array with values the same as described
%   in the documentation for WRITETOLOG below.

if ( nargin < 2 )
    quietMode = false;
end

% Store the current directory: (We will return here when we're done.)
currDir = pwd();

% Switch to the test directory and get list of all the test files:
cd(testDir);
testFiles = dir('*.m');
testFiles = {testFiles.name};

% For making the output string align nicely:
maxLength = max(cellfun(@length, testFiles));

% Allocate pass and timing variables:
numFiles  = numel(testFiles);
durations = zeros(numFiles, 1);
errorMessages = {'FAILED', 'CRASHED'};

% Attempt to run all of the tests:
try % Note, we try-catch as we've CD'd and really don't want to end up elsewhere
    
    % Loop over the test files:
    for k = 1:numFiles
        % Next file to test: (.m extension is removed).
        testFile = testFiles{k}(1:end-2);
        
        if ( quietMode )
            % --quiet mode
            durations(k) = runTest(testFile);
            if ( durations(k) < 0 )
                printTestInfo(testDir, testFile, k, maxLength);
                message = errorMessages{-durations(k)};
                fprintf([message '\n']);
            end
            
        else
            % --verbose mode
            printTestInfo(testDir, testFile, k, maxLength);
            durations(k) = runTest(testFile);
            if ( durations(k) > 0 )
                % Success message.
                message = sprintf('passed in %.4fs', durations(k));
            else
                % Error flags are negative indices.
                message = errorMessages{-durations(k)};
            end
            fprintf([message '\n']);

        end

    end
    
catch ME
    % We failed. Return to the starting directory and rethrow the error:
    cd(currDir)
    rethrow(ME)
    
end

% Yay!
if ( all(durations > 0) )
    fprintf('All %s tests passed in %4.4fs.\n', testDir, sum(durations));
elseif ( ~quietMode )
    % Note. We don't show this in quiet mode as it's already clear.
    fprintf('%d failed test(s) in %s directory.\n', sum(durations < 0), testDir);
end

% Restore the current working directory and return:
cd(currDir);

% Dump all the test data into a cell array to pass back to the CHEBTEST
% function. This data is what is written to a .CSV file if logging is on.
testResults = cell(numFiles,3);
for k = 1:numFiles
    testResults{k,1} = testDir;               % directory name
    testResults{k,2} = testFiles{k}(1:end-2); % test file name
    testResults{k,3} = durations(k);          % duration / error flag
end

end


function printTestInfo(testDir, testFile, k, maxLength)
%PRINTTESTINFO   Pretty print test info for a given test file.
%   PRINTTESTINFO(TESTDIR, TESTFILE, K) will print the following to screen:
%     Test #k: TESTDIR/TESTFILE.m ...
%   Note that it will not linebreak after this.
%
%   PRINTTESTINFO(TESTDIR, TESTFILE, K, MAXLENGTH) is similar, but will print
%   extra whitespace after the ellipsis so that the following text for each
%   TESTFILE in TESTDIR is aligned.
%
%   Furthermore, PRINTTESTINFO tests to see if the Matlab desktop is running and
%   will support HTML. If it is, the displayed the displayed TESTDIR/TESTFILE.m
%   is hyperlinked to open TESTFILE.m in the Matlab editor.

link = printTestFilename(testDir, testFile);           % Either html or plain.
ws = repmat(' ', 1, maxLength - length(testFile) - 1); % Whitespace.
fprintf('  Test #%3.3d: %s ...%s', k, link, ws);       % Print to screen.

end

function link = printTestFilename(testDir, testFile, loc)
%PRINTTESTFILENAME   Pretty print a filename.
%   STR = PRINTTESTFILENAME(TESTDIR, TESTFILE, LOC) tests to see if the Matlab
%   desktop is running and will support HTML. If it is, the displayed the
%   displayed TESTDIR/TESTFILE.m is hyperlinked to open LOC/TESTFILE.m in the
%   Matlab editor. If LOC is not passed, PRINTTESTFILENAME will attempt to find
%   the location of TESTFILE.M via fileparts(which(TESTFILE)).

if ( usejava('jvm') && usejava('desktop') )
    % Use HTML links if Java is enabled. 
    
    if ( nargin < 3 )
        loc = fileparts(which(testFile));
    end
    url = fullfile(loc, [testFile '.m']);
    
    link = ['<a href = "matlab: edit ''', url, '''">' testDir '/', ...
        testFile '.m</a>'];  
    
else
    
    % Otherwise, use plaintext.
    link = fullfile(testDir, [testFile '.m']);
    
end

end


function duration = runTest(testFile)
%RUNTEST Runs the given test file.
%   DURATION = RUNTEST(TESTFILE) executes the file TESTFILE.
%   TESTFILE should return a vector or matrix of logical values. 
%
%   If each of these are logical true (or 1) then the test passes and RUNTEST
%   returns DURATION as a double > 0 corresponding to the time it took the test
%   to execute in seconds.
%
%   If any of the entries returned by TESTFILE are logcial false (or 0), then
%   DURATION = -1.
%
%   If executing TESTFILE crashes, this is caught in a try-catch statement, and
%   RUNTTEST returns DURATION = -2.

% Close any open windows:
close all

% Attempt to run the test:
try
    tstart = tic();
    pass = feval(testFile);
    duration = toc(tstart);

    pass = all(pass(:));
    if ( ~pass )
        % The test failed, so return FAILED flag.
        duration = -1;
    end

catch ME %#ok<NASGU>
    % We crashed. This is bad. Return CRASHED flag.
    duration = -2;
    
    % But we _don't_ want to throw an error.
    %rethrow(ME)
end

% Close any windows the test may have left open:
close all

end
