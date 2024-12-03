function out = kdmroot()
%KDMROOT   Root directory of Kerenl DM Suite installation.
%   S = KDMROOT() returns a string that is the name of the directory where
%   Kernel DM Suite is installed.

out = fileparts(which('kdmroot'));

end