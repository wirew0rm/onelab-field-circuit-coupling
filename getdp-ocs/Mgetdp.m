function [a,b,c] = Mgetdp (string,parameters,parameternames,extvar,intvar,t)

  verbosity = [0;0];
	diffdelta = 1e-5;       % relative finite difference perturbation, disabled if zero
  
  % check dimensions
  if (length(extvar)~=2)
    error('only models with a single inductance supported.')
  end
  
  % check file name
  if ~exist([string '/' string '.pro'],'file')
    error('file "%s" dos not exist.',string);
  end

  % init variables
  if isempty(intvar)
    intvar = zeros(2,1);
  end
  u   = extvar;
  phi = intvar(1);
  jL  = intvar(2);
  
  % set a minimal excitation to avoid numerical difficulties
  % inductance L will be constant for small values jL and thus this 
  % workaround will not affect the solution
  if jL<1e-10
    jgetdp = 1e-10;
  else
    jgetdp = jL;
  end
  
  % start computation in GetDP
	% set coil excitation
	Iparam = ol_getParameters('Input/4Coil Parameters/0Current (rms) [A]');
	Iparam{1}.value = jgetdp;
	ol_setParameter(Iparam{1});
	% run GetDP
	ol_runBlockingSubClient('incuctance',['getdp ' string '/' string '.pro -pre -cal -pos -msh ' string '/' string '.msh -gmshread ' string '/res/a.pos']);
	% get inductance
	Lparam = ol_getParameters('Output/50Inductance from Flux [H]');
  L = Lparam{1,1}.value;
	
	Initparam = ol_getParameters('Input/42InitSolutionFromPrevious');
	Initparam{1}.value = 1;
	ol_setParameter(Initparam{1});

	dLdjL = 0;
	if diffdelta > 0
		% set coil exitation
		Iparam = ol_getParameters('Input/4Coil Parameters/0Current (rms) [A]');
		Iparam{1}.value = jgetdp + jgetdp*diffdelta;
		ol_setParameter(Iparam{1});
		% run getdp
		ol_runBlockingSubClient('incuctance',['getdp ' string '/' string '.pro -pre -cal -pos -msh ' string '/' string '.msh -gmshread ' string '/res/a.pos']);
		% get inductance
		Lparam = ol_getParameters('Output/50Inductance from Flux [H]');
		Ld = (Lparam{1,1}.value*(1+diffdelta) - L)/diffdelta;
	else
		% simply use L
		Ld = L;
	end
  
  % debug output
	printf(['Mgetdp: L = ' num2str(L,'%1.5e') 'H ,Ld = ' num2str(Ld,'%1.5e') 'H, I = ' num2str(jgetdp,'%1.5e') 'A at t=' num2str(t,'%1.1e') 's.\n']);
  ol_sendInfo(['Mgetdp: L = ' num2str(L,'%1.5e') 'H ,Ld = ' num2str(Ld,'%1.5e') 'H, I = ' num2str(jgetdp,'%1.5e') 'A at t=' num2str(t,'%1.1e') 's.\n']);
  
  a = [0 0 0 0; 
       0 0 0 0; 
       0 0 1 0; 
       0 0 0 0];

  % b should not contain L but dL/djL 
  b = [0 0 0  1; 
       0 0 0 -1; 
      -1 1 0  0; 
       0 0 1 -Ld];

  % c should indeed contain phi-L*jL 
  d = [0 0 0  1; 
       0 0 0 -1; 
      -1 1 0  0; 
       0 0 1 -L];

  c = d*[u;phi;jL];
  
end %function
