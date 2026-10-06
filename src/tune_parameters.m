function [bestV,bestA,logTable] = tune_parameters(prob,cfg)
%TUNE_PARAMETERS Empirical tuning; every evaluation logs parameters and seed.
% Fixed-budget tuning disables tolerance stopping. Stability is empirical,
% not a convergence guarantee. Failure raises an error and saves the log.
records=struct('method',{},'stage',{},'round',{},'candidate',{},...
    'budget',{},'seed',{},'cAlpha',{},'eta',{},'cTheta',{},...
    'cGamma',{},'cBeta',{},'theta2',{},'score',{},'finalRSE',{},...
    'iterations',{},'stopReason',{},'failed',{});
v=arrayfun(@(x) struct('cAlpha',x),cfg.tuning.vrrkCAlpha,...
    'UniformOutput',false);
bestV=race(v,"VRRK","stepsize");
a={};
for scale=cfg.tuning.avStepScale
    for eta=cfg.tuning.eta
        a{end+1}=av(bestV.cAlpha*scale,eta,1,1,0.5); %#ok<AGROW>
    end
end
a1=race(a,"AVRRK","stepsize_eta");
a={};
for ct=cfg.tuning.cTheta
    for cg=cfg.tuning.cGamma
        a{end+1}=av(a1.cAlpha,a1.eta,ct,cg,0.5); %#ok<AGROW>
    end
end
a2=race(a,"AVRRK","momentum");
a=cell(numel(cfg.tuning.theta2),1);
for j=1:numel(a)
    a{j}=av(a2.cAlpha,a2.eta,a2.cTheta,a2.cGamma,cfg.tuning.theta2(j));
end
bestA=race(a,"AVRRK","theta2");
% Both methods must pass long-horizon validation.
bestV=stability(bestV,"VRRK");
bestA=stability(bestA,"AVRRK");
logTable=struct2table(records);

    function p=av(ca,eta,ct,cg,t2)
        p=struct('mode',"empirical",'cAlpha',ca,'eta',eta,...
            'cTheta',ct,'cGamma',cg,'cBeta',1,'theta2',t2);
    end

    function h=evaluate(p,method,stage,round,ci,budget,seed)
        c=cfg;
        c.tol=0;
        % Demo iteration limits must not truncate tuning/validation budgets.
        c.maxIterations=Inf;
        c.maxPasses=budget;
        rng(seed,'twister');
        if method=="VRRK"
            h=run_vrrk_as(prob,p,c);
        else
            h=run_avrrk_as(prob,p,c);
        end
        row=struct('method',method,'stage',stage,'round',round,...
            'candidate',ci,'budget',budget,'seed',seed,...
            'cAlpha',p.cAlpha,'eta',NaN,'cTheta',NaN,'cGamma',NaN,...
            'cBeta',NaN,'theta2',NaN,'score',history_score(h),...
            'finalRSE',h.rse(end),'iterations',h.iter(end),...
            'stopReason',h.stopReason,'failed',h.failed);
        fields={'eta','cTheta','cGamma','cBeta','theta2'};
        for fieldIndex=1:numel(fields)
            name=fields{fieldIndex};
            if isfield(p,name)
                row.(name)=p.(name);
            end
        end
        records(end+1)=row;
    end

    function best=race(candidates,method,stage)
        active=1:numel(candidates);
        for round=1:numel(cfg.tuning.budgets)
            budget=cfg.tuning.budgets(round);
            scores=inf(numel(active),1);
            for u=1:numel(active)
                ci=active(u);
                vals=inf(2,1);
                for sj=1:2
                    h=evaluate(candidates{ci},method,stage,round,ci,...
                        budget,cfg.seed+10000*round+sj);
                    if ~h.failed && h.rse(end)<h.rse(1)
                        vals(sj)=history_score(h);
                    end
                end
                % A candidate must remain finite on both validation streams.
                if all(isfinite(vals))
                    scores(u)=median(vals);
                end
            end
            [scores,order]=sort(scores);
            order=order(isfinite(scores));
            if isempty(order)
                % Revisit ALL original candidates at the longest budget.
                % Short-budget elimination can discard stable slow starters.
                best=recover(candidates,method,stage);
                return;
            end
            if round<numel(cfg.tuning.budgets)
                keep=min(numel(order),max(1,ceil(numel(active)*cfg.tuning.keepFraction)));
                order=order(1:keep);
            end
            active=active(order);
        end
        best=candidates{active(1)};
    end

    function best=recover(candidates,method,stage)
        % Bounded backtracking; never accept a failed or non-improving run.
        % Scale 1 first retries candidates eliminated by successive halving.
        scales=2.^(-(0:8));
        budget=max(cfg.tuning.budgets);
        fprintf('%s: recovering stage %s at budget %g.\n',method,stage,budget);
        for attempt=1:numel(scales)
            scores=inf(numel(candidates),1);
            for ci=1:numel(candidates)
                p=candidates{ci};
                p.cAlpha=p.cAlpha*scales(attempt);
                vals=inf(2,1);
                for sj=1:2
                    h=evaluate(p,method,stage+"_recovery",attempt,ci,...
                        budget,cfg.seed+10000*numel(cfg.tuning.budgets)+sj);
                    if ~h.failed && h.rse(end)<h.rse(1)
                        vals(sj)=history_score(h);
                    end
                end
                if all(isfinite(vals))
                    scores(ci)=median(vals);
                end
            end
            [score,ci]=min(scores);
            if isfinite(score)
                best=candidates{ci};
                best.cAlpha=best.cAlpha*scales(attempt);
                fprintf('%s: recovered %s with cAlpha=%g.\n',method,stage,best.cAlpha);
                return;
            end
        end
        fail(method+" has no improving candidate after bounded recovery in "+stage);
    end

    function best=stability(initial,method)
        scales=2.^(-(0:8));
        scores=inf(size(scales));
        candidates=cell(size(scales));
        for ci=1:numel(scales)
            p=initial;
            p.cAlpha=p.cAlpha*scales(ci);
            candidates{ci}=p;
            terminal=inf(2,1);
            for sj=1:2
                h=evaluate(p,method,"stability",0,ci,...
                    max([cfg.tuning.stabilityPasses cfg.tuning.budgets]),cfg.seed+777+sj);
                if isfinite(history_score(h)) && h.rse(end)<h.rse(1)
                    terminal(sj)=h.rse(end);
                end
            end
            if all(isfinite(terminal))
                scores(ci)=median(terminal);
            end
        end
        [score,ci]=min(scores);
        if ~isfinite(score)
            fail(method+" has no candidate passing the stability check");
        end
        best=candidates{ci};
    end

    function fail(message)
        failureLog=struct2table(records);
        failureFile=[tempname '.mat'];
        % Save the training problem as well, so a failure is reproducible.
        save(failureFile,'failureLog','cfg','prob');
        error('VRRK:TuningFailed','%s. Evaluation log saved to %s.',message,failureFile);
    end
end
