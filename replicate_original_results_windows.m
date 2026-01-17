% replicate_original_results.m  (Windows-path version)
% Deterministic runner for the author's MATLAB code on your Windows machine.
% - cd to user-specified code directory
% - rng(123)
% - diary to replication_outputs\logs
% - run('main.m')
% - save workspace to replication_outputs\run_workspace.mat

function replicate_original_results()
    try
        %% Change to your Windows code directory
        cd('C:\Users\arpit\OneDrive - University of Pittsburgh\Pitt PhD Study Material\11 RESEARCH\Second Year\Marla Solving in Levels\Code\Code for JME');
        
        %% Deterministic seed
        rng(123);
        
        %% Ensure output folders exist
        outdir = 'C:\Users\arpit\OneDrive - University of Pittsburgh\Pitt PhD Study Material\11 RESEARCH\Second Year\Marla Solving in Levels\replication_outputs';
        logdir = 'C:\Users\arpit\OneDrive - University of Pittsburgh\Pitt PhD Study Material\11 RESEARCH\Second Year\Marla Solving in Levels\replication_outputs\logs';
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        if ~exist(logdir, 'dir'), mkdir(logdir); end
        
        %% Start diary
        diary(fullfile(logdir, 'run_20250921_215113.txt'));
        diary on;
        fprintf('\n=== Replication Run Started: %s ===\n', datestr(now));
        fprintf('CWD: %s\n', pwd);
        fprintf('Seed: 123\n');
        
        %% (Optional) Centralize solver options
        % opts = optimoptions('fsolve','Display','iter','TolFun',1e-8,'TolX',1e-8,'MaxIterations',5000);
        % Pass opts into each fsolve call in subroutines if feasible.
        
        %% Run the author's main script
        run('main.m');
        
        %% Save workspace snapshot
        save(fullfile(outdir, 'run_workspace.mat'));
        
        %% Basic budget check (only if variables are present)
        try
            if exist('Pc','var') && exist('C','var') && exist('Px','var') && exist('X','var') ...
               && exist('NX','var') && exist('w','var') && exist('L','var') ...
               && exist('r','var') && exist('K','var')
                bc_resid = max(abs(Pc.*C + Px.*X + NX - (w.*L + r.*K)), [], 'all');
                fprintf('\nMax budget-constraint residual: %g\n', bc_resid);
            else
                fprintf('\n[Note] Budget variables not all found; skipping BC residual check.\n');
            end
        catch ME
            fprintf('\n[Warning] Post-run BC check failed: %s\n', ME.message);
        end
        
        diary off;
        fprintf('\n=== Replication Run Finished ===\n');
    catch ME
        try, diary off; end
        rethrow(ME);
    end
end
