#!/usr/bin/octave-cli -q
%# import packages
pkg load ocs       % circuit simulation
pkg load fides     % additions to ocs for coupling
pkg load onelab    % octave interface to onelab

% parameter parsing
for i = [1:nargin]
  if (i+2 > nargin)
    error("invoke with -onelab name server_addr");
  elseif (strcmp(argv(){i},"-onelab"))
    name = argv(){i+1};
    addr = argv(){i+2};
    break;
  endif
endfor

% Start the client
ol_client(name, addr);

% check which Action is requested by the UI and run the according script.
action = ol_getParameters([name '/Action']);
if     (!isempty(action) && strcmp(action{1}.value, 'compute'))
	ocs_onelab_compute;
elseif (!isempty(action) && strcmp(action{1}.value, 'initialize'))
	ocs_onelab_init;
elseif (!isempty(action) && strcmp(action{1}.value, 'check'))
	ocs_onelab_init;
elseif (!isempty(action))
	ol_sendInfo(['Nothing to do: ' action{1}.value]);
else 
	ol_sendInfo('Nothing to do');
endif

% disconnect from onelab server and exit
ol_disconnect();
exit;
