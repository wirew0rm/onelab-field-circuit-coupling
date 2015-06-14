% get Parameters from onelab
tmp= ol_getParameters('rl/e1');
e1 = tmp{1};
tmp= ol_getParameters('rl/e2'); 
e2 = tmp{1};
tmp= ol_getParameters('rl/tmax');
tmax = tmp{1};
tmp= ol_getParameters('rl/t');
t = tmp{1};

%initialize ocs
cirstruct = prs_iff ('rl/rl');                                                % init circuit model
T         = [t.value tmax.value];                                             % time interval
hinit     = 1e-8;                                                             % initial step size
tout      = T(1):hinit:T(end);                                                % interval: fixed steps
x0        = [0;0;0;0;0];                                                      % initial value
dmp       = 1;                                                              % damping parameter for newton
ntol      = 1e-4;                                                             % tolerance for newton
maxit     = 40;                                                               % maximal iterations of newton
pltvars   = {'e_1','e_2'};                                                    % plot variables
idxvars   = utl_index_by_name(cirstruct,pltvars);                             % indices of pltvars
verbosity = [1 0];                                                            % verbosity level (1=on, 0=off)
showplot  = 0;                                                                % set >0 for plot of the result

%# time integration
ol_sendInfo(['Time integration ' num2str(length(tout)) ' timesteps']);
printf(['Time integration ' num2str(length(tout)) ' timesteps\n']);
tic;                                                                          % start stop watch
[xout] = tst_backward_euler_noguess(cirstruct,x0,tout,...                     % fixed time stepping with backward euler
         ntol,maxit,pltvars,verbosity);
cpu = toc                                                                     % remember the used cpu time

% post results to onelab
e1.choices = num2cell(xout(idxvars(1),:));
e1.value = xout(idxvars(1),end);
ol_setParameter(e1);
e2.choices = num2cell(xout(idxvars(2),:));
e2.value = xout(idxvars(2),end);
ol_setParameter(e2);
t.choices = num2cell(tout(:));
t.value = tout(end);
ol_setParameter(e3);

%# plot solution
if showplot
  figure;
  plot(tout,xout(idxvars,:),'-');
  title('A simple RL circuit');
  legend(pltvars);
end
