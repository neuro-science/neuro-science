function evt = cmp16_rawmeg_evt(evr, sampExt)
% % % 13/06/16	updated by wp to add more compatibility: 
% % % 25/07/14	written by wp: 

	%% check input
	if nargin < 2
		sampExt = [180, 12000];
	end
	
	%% get all events
	evr = struct2cell(evr);	% get original data events
	nTriggers = 0;		%initial number of triggers and responses
	nResponses = 0;
	trgType = zeros(2, 1) + nan;
	trgPointsOrigin = zeros(2, 1) + nan;
	rspType = zeros(2, 1) + nan;
	rspPointsOrigin = zeros(2, 1) + nan;
	for iEvt = 1 : size(evr, 2)
		switch evr{1, iEvt}		% sort events
			case 'UPPT001'		%trigger of stimulus
				nTriggers = nTriggers + 1;
				trgType(nTriggers) = evr{3, iEvt};
				trgPointsOrigin(nTriggers) = evr{2, iEvt};
			case 'UPPT002'		%response
				nResponses = nResponses + 1;
				rspType(nResponses) = evr{3, iEvt};
				rspPointsOrigin(nResponses) = evr{2, iEvt};
			otherwise;
		end
	end
	% % % if there is no triggers at all
	if nTriggers < 1
		evt = [];
		return;
	end
	
	rspTypeKey = log2(rspType - 224) * 2 - 5;	%key to key id
	rspTypeKey(abs(rspTypeKey) > 1) = 0;
	
	%% get all events
	evt.nEvents = nTriggers;	%initial event data
	evt.StiType = zeros(evt.nEvents, 1) + nan;
	evt.TrigSample = zeros(evt.nEvents, 1) + nan;
	evt.RspSample = zeros(evt.nEvents, 1) + nan;
	evt.RspType = zeros(evt.nEvents, 1) + nan;
	evt.RspSamplePools = cell(nTriggers, 1);	%size  nTriggers
	evt.RspTypePools = cell(nTriggers, 1);	%size  nTriggers
	for iTrg = 1 : evt.nEvents
		if iTrg < nTriggers
			t0 = trgPointsOrigin(iTrg + 1);
		else
			t0 = trgPointsOrigin(iTrg) + sampExt(2);
		end
		evt.StiType(iTrg) = trgType(iTrg);	%set value
		evt.TrigSample(iTrg) = trgPointsOrigin(iTrg);	%set value
		idx = (rspPointsOrigin > trgPointsOrigin(iTrg)) & (rspPointsOrigin < t0);
		evt.RspSamplePools{iTrg} = (rspPointsOrigin(idx) - evt.TrigSample(iTrg));
		evt.RspTypePools{iTrg} = rspTypeKey(idx);
		n = length(evt.RspSamplePools{iTrg});
		for k = 1 : n
			if evt.RspSamplePools{iTrg}(k) > sampExt(1);
				evt.RspSample(iTrg) = evt.RspSamplePools{iTrg}(k);
				evt.RspType(iTrg) = evt.RspTypePools{iTrg}(k);
				break;
			end
		end
	end
end %end of function

